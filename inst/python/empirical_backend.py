#!/usr/bin/env python3
"""Portable paper-specific empirical workflows; never downloads market data."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"

try:
    import numpy as np
    import pandas as pd
    import scipy
    from scipy import stats
    from statsmodels.stats.diagnostic import acorr_ljungbox, het_arch
    from statsmodels.stats.multitest import multipletests
    from threadpoolctl import threadpool_limits
except ImportError as exc:
    raise SystemExit("Empirical backend dependency missing: " + str(exc) +
                     ". Select a Python with numpy, pandas, scipy, statsmodels and threadpoolctl installed.")

import empirical_core as core
import empirical_data as data
import break_core as breaks_core
import quic_solver


def serial(value):
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def save_json(value, path):
    Path(path).write_text(json.dumps(serial(value), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def validate(y, f):
    y, f = np.asarray(y, dtype=float), np.asarray(f, dtype=float)
    if y.ndim != 2 or f.ndim != 2 or len(y) != len(f):
        raise ValueError("Y and F must be numeric matrices with the same number of rows")
    if y.shape[1] < 40 or len(y) <= f.shape[1] + 5:
        raise ValueError("The full grid needs N >= 40 and residual degrees of freedom > 4")
    if not np.isfinite(y).all() or not np.isfinite(f).all():
        raise ValueError("Y and F must be finite, complete matrices")
    design = np.column_stack([np.ones(len(y)), f])
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("F plus an intercept must have full column rank; do not supply an intercept column")
    return core.fit(y, f)


def alpha_table(a, assets):
    p = 2 * stats.t.sf(np.abs(a["t"]), a["v"])
    return pd.DataFrame(dict(asset=assets, alpha_monthly=a["coef"][0],
                             alpha_bp_monthly=10000*a["coef"][0],
                             standard_error_bp=10000*np.sqrt(a["s2"]/a["q"]),
                             standardized_alpha=a["t"], individual_p=p,
                             BH_q=multipletests(p, method="fdr_bh")[1]))


def empirical_test(y, f, B, seed, include_maxk, B_maxk, output_dir, assets=None,
                   setting="custom", maxk_seed=None, consume_l_before_maxk=False):
    """Same Rademacher stream and batch-128 arithmetic as the final oldstream run."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    quic_solver.configure(output_dir / "compiled")
    a = validate(y, f)
    grid = core.grid(a["N"])
    orders = np.array(sorted({1, *grid}), dtype=int)
    position = np.searchsorted(orders, grid)
    observed = core.ordered_sums(a["t"]**2, orders)
    boot, arithmetic = core.original_stream_bootstrap(a, orders, B, seed)
    result = core.evaluate(observed[position], boot[:, position], a["N"], a["T"])
    py_stat, py_p, _, _ = core.py_test(a)
    mx = float(np.max(a["t"]**2))
    gm = mx - 2*np.log(a["N"]) + np.log(np.log(a["N"]))
    max_p = float(-np.expm1(-np.exp(-gm/2)/np.sqrt(np.pi)))
    rows = [dict(setting=setting, method="AdaL", statistic=result["score"], p_value=result["p"]),
            dict(setting=setting, method="T1", statistic=observed[0], p_value=(1+np.sum(boot[:, 0] >= observed[0]))/(B+1)),
            dict(setting=setting, method="PY", statistic=py_stat, p_value=py_p),
            dict(setting=setting, method="MAX", statistic=mx, p_value=max_p),
            dict(setting=setting, method="COM", statistic=np.nan, p_value=min(1., 2*min(py_p, max_p)))]
    maxk_meta = None
    if include_maxk:
        if B_maxk < 2:
            raise ValueError("max-k requires B_maxk >= 2 to estimate bootstrap standard deviations")
        rng = np.random.default_rng(seed if maxk_seed is None else maxk_seed)
        if consume_l_before_maxk:
            # Final China comparator follows its 9999-draw bimarket L stream,
            # distinct from the 99999-draw AdaL stream retained in the final paper.
            for start in range(0, B_maxk, 128):
                rng.choice([-1., 1.], size=(min(128, B_maxk-start), a["T"]))
        xr, maxk_meta, _, _ = core.xia_test(a, B_maxk, rng)
        for k in (5, 10, 30):
            rows.append(dict(setting=setting, method="max-k (K=%d)" % k,
                             statistic=xr["maxk_K%d_stat" % k], p_value=xr["maxk_K%d_p" % k]))
    table = pd.DataFrame(rows)
    table["B"] = [B if method in ("AdaL", "T1") else
                  B_maxk if method.startswith("max-k") else np.nan for method in table.method]
    table["N"], table["T"] = a["N"], a["T"]
    table.to_csv(output_dir / "results.csv", index=False)
    components = pd.DataFrame(dict(k=grid, statistic=observed[position],
                                   exceedances=result["exceedance"], p_value=result["component_p"],
                                   clipped_p=result["clipped_p"], weight=1/len(grid)))
    components.to_csv(output_dir / "components.csv", index=False)
    if assets is None:
        assets = [str(i+1) for i in range(a["N"])]
    alphas = alpha_table(a, assets)
    alphas.to_csv(output_dir / "alphas.csv", index=False)
    save_json(dict(setting=setting, seed=seed, rng="numpy.random.default_rng / PCG64", B=B,
                   batch_size=128, N=a["N"], T=a["T"], residual_df=a["v"],
                   component_calibration="add-one empirical upper tail for every full-grid component",
                   denominator="refitted for every wild-bootstrap draw", endpoint_epsilon=result["epsilon"],
                   COM="min(1, 2 * min(PY_p, MAX_p))", maxk_B=B_maxk if include_maxk else None,
                   maxk_seed=maxk_seed, maxk_follows_L_stream=consume_l_before_maxk,
                   arithmetic=arithmetic, maxk_solver=maxk_meta), output_dir / "metadata.json")
    return table, components, alphas


def load_settings(data_dir, output_dir, requested):
    data.configure(data_dir, Path(output_dir)/"data_audit")
    settings = []
    if "china" in requested:
        raw = data.get_data()
        y, f, attr, _ = data.panel(raw, dict(start="2021-01", end="2025-12", pool="top1000"))
        settings.append(dict(key="china", y=y.to_numpy(), f=f.to_numpy(), assets=y.columns,
                             months=y.index, base_f=f.to_numpy(), break_date=None, attr=attr,
                             normalized=(1+raw["returns"].reindex(index=y.index, columns=y.columns)).cumprod(),
                             B=99999, seed=core.seed_from_legacy(20260903, "2021-01_2025-12_top1000alpha1"),
                             maxk_seed=core.seed_from_legacy(20260906, "CN_main_2021-01_2025-12_top1000alpha1"),
                             consume_l_before_maxk=True))
    if any(s in requested for s in ("us_ff3", "us_ff5")):
        raw = data.read_us_monthly()
        y, f3, attr, assets = data.us_panel(raw, "2012-01", "2019-12", "continuous")
        for key, ff5, break_date, seed_key, maxk_key in (
                ("us_ff3", False, "2017-12", "US_piecewise_loading_2012_2019alpha1", "US_piecewise_loading_2012_2019_maxk"),
                ("us_ff5", True, "2015-07", "US_FF5_piecewise_L_AdaL", "US_FF5_piecewise_maxk")):
            if key not in requested:
                continue
            yy, f = y.copy(), f3.copy()
            if ff5:
                ff = data.french(Path(data_dir)/"America/F-F_Research_Data_5_Factors_2x3.csv", five=True).reindex(y.index)
                if ff.isna().any().any():
                    raise ValueError("Missing monthly observations in FF5 factors")
                yy = yy.add(raw["factors"]["RF"].reindex(y.index), axis=0).sub(ff["RF"], axis=0)
                f = ff[["MKT", "SMB", "HML", "RMW", "CMA"]]
            base_f = f.to_numpy().copy()
            post = (f.index >= pd.Period(break_date, "M")).astype(float)
            f = f.copy()
            for column in list(f.columns):
                f[column+"_POST"] = f[column]*post
            settings.append(dict(key=key, y=yy.to_numpy(), f=f.to_numpy(), base_f=base_f,
                                 months=y.index, assets=assets, break_date=break_date, attr=attr,
                                 normalized=raw["prices"].reindex(index=y.index, columns=assets),
                                 B=19999, seed=core.seed_from_legacy(20260906, seed_key),
                                 maxk_seed=core.seed_from_legacy(20260906, maxk_key), consume_l_before_maxk=False))
    for item in settings:
        n_expected, t_expected = (980, 60) if item["key"] == "china" else (288, 96)
        if item["y"].shape != (t_expected, n_expected):
            raise ValueError("%s paper sample should have shape (%d, %d); obtained %s. Check input vintages and sample-selection fields. Use lfpm_empirical_test for a different sample."
                             % (item["key"], t_expected, n_expected, item["y"].shape))
        validate(item["y"], item["f"])
    return settings


def penalty_sensitivity(items, output_dir, B, seed_base=20260906):
    """Recompute legacy constant-loading penalty sensitivities with explicit seeds.

    The historical summary did not archive these RNG seeds. This interface
    makes the new deterministic convention explicit instead of inventing parity.
    """
    if B < 2:
        raise ValueError("Sensitivity max-k requires at least two draws")
    destination = output_dir/"sensitivity"
    destination.mkdir(parents=True, exist_ok=True)
    quic_solver.configure(destination/"compiled")
    rows = []
    for item in items:
        if item["key"] == "us_ff5":
            continue
        a = validate(item["y"], item["base_f"])
        panel = "CN_main_2021-01_2025-12_top1000" if item["key"]=="china" else "US_roll_2012-01_2019-12_T96"
        for scale in (.75, 1., 1.25):
            key = panel+"rho"+str(scale)
            seed = core.seed_from_legacy(seed_base, key)
            result, meta, _, _ = core.xia_test(a, B, np.random.default_rng(seed), rho_scale=scale)
            for k in (5, 10, 30):
                rows.append(dict(setting=item["key"], model="constant_loadings", rho_scale=scale,
                    candidate_orders=k, statistic=result["maxk_K%d_stat" % k],
                    p_value=result["maxk_K%d_p" % k], B=B, seed=seed, seed_key=key,
                    historical_rng_parity_verified=False, **meta))
    pd.DataFrame(rows).to_csv(destination/"maxk_penalty_sensitivity.csv", index=False)
    save_json(dict(seed_base=seed_base,
        status="recomputed_with_explicit_seed_convention",
        limitation="Original bimarket penalty-sensitivity RNG seeds were not archived in the available source. The reported 0.0264 endpoint has not been verified by this new run.",
        grids=[.75,1.,1.25], B=B), destination/"provenance.json")


def residual_diagnostics(item, out):
    a = validate(item["y"], item["f"])
    rows = []
    for j, asset in enumerate(item["assets"]):
        lb = acorr_ljungbox(a["e"][:, j], lags=[6, 12], boxpierce=True, return_df=True)
        rows.append(dict(asset=asset, BP6=lb.bp_pvalue.iloc[0], BP12=lb.bp_pvalue.iloc[1],
                         LB12=lb.lb_pvalue.iloc[1], ARCH6=het_arch(a["e"][:, j], nlags=6, ddof=a["X"].shape[1])[1]))
    frame = pd.DataFrame(rows)
    summary = []
    for name in ("BP6", "BP12", "LB12", "ARCH6"):
        q = multipletests(frame[name], method="fdr_bh")[1]
        summary.append(dict(setting=item["key"], diagnostic=name, raw_count=int((frame[name]<.05).sum()), BH_count=int((q<=.05).sum())))
        frame[name+"_BH_q"] = q
    frame.to_csv(out/"residual_diagnostics.csv", index=False)
    pd.DataFrame(summary).to_csv(out/"diagnostic_summary.csv", index=False)


def model_diagnostics(item, out, B=1999):
    """Constant-loading specification scan; descriptive residual gap ratio."""
    fitted = validate(item["y"], item["f"])
    e, half = fitted["e"], fitted["T"]//2
    transformed = e[:half] + 1j*e[half:2*half]
    eigenvalues = np.linalg.eigvalsh(transformed@transformed.conj().T/half)[::-1].real
    gaps = eigenvalues[:7]-eigenvalues[1:8]
    gap_ratio = float(np.max(gaps[:5]/np.maximum(gaps[1:6], 1e-14)))
    # This scan is on the constant-loading model, including for U.S. FF3.
    a = validate(item["y"], item["base_f"])
    T, N = a["T"], a["N"]
    X, Q, es = a["X"], a["Q"], a["e"]/np.sqrt(a["s2"])
    candidates = np.arange(max(12, int(.2*T)), min(T-12, int(.8*T))+1)
    operators = []
    for t0 in candidates:
        h = X[:, 1:]*(np.arange(T)>=t0)[:, None]
        h = h-Q@(Q.T@h)
        operators.append(np.linalg.qr(h)[0].T)
    flat = np.asarray(operators).reshape(-1, T)
    p = item["base_f"].shape[1]
    squared = ((flat@es).reshape(len(candidates), p, N)**2).sum(axis=1)
    observed = np.array([squared.mean(axis=1).max(), squared.max()])
    seed_key = "US_FF5_loading_break_scan" if item["key"]=="us_ff5" else "slope"+str(T)+str(N)
    seed = core.seed_from_legacy(20260906, seed_key)
    rng = np.random.default_rng(seed)
    boot = []
    for _ in range(B):
        signs = rng.choice([-1., 1.], T)
        score = (flat@(es*signs[:, None])).reshape(len(candidates), p, N)
        score2 = (score**2).sum(axis=1)
        boot.append([score2.mean(axis=1).max(), score2.max()])
    probability = (1+(np.asarray(boot)>=observed).sum(axis=0))/(B+1)
    selected = int(candidates[squared.mean(axis=1).argmax()])
    pd.DataFrame([dict(setting=item["key"], static_complex_gap_ratio=gap_ratio,
        loading_average_p=probability[0], loading_maximum_p=probability[1],
        loading_bonferroni_p=min(1., 2*probability.min()), loading_B=B,
        loading_seed=seed, loading_average_break=str(item["months"][selected]),
        gap_ratio_is_descriptive=True)]).to_csv(out/"model_diagnostics.csv", index=False)


def fitted_null(items, output_dir, repetitions, inner_B):
    """Recompute both denominator choices on the original paired seed streams."""
    rows, comp = [], []
    for item in items:
        if item["key"] == "us_ff5":
            continue
        base = validate(item["y"], item["f"])
        order = core.grid(base["N"])
        fitted = base["X"][:, 1:] @ base["coef"][1:]
        for rep in range(repetitions):
            if item["key"] == "china":
                outer = core.seed_from_legacy(20260903, "outer%d" % rep)
                inner = core.seed_from_legacy(20260903, "sim%dalpha1" % rep)
            else:
                outer, inner = 20260928730100+rep, 20260928730200+rep
            signs = np.random.default_rng(outer).choice([-1., 1.], size=base["T"])
            a = core.fit(fitted+base["e"]*signs[:, None], item["f"])
            obs = core.ordered_sums(a["t"]**2, order)
            frozen = core.bootstrap(a, order, inner_B, inner)
            refitted = np.empty_like(frozen)
            rng = np.random.default_rng(inner)
            ss = np.sum(a["e"]**2, axis=0)
            for start in range(0, inner_B, 256):
                count = min(256, inner_B-start)
                xi = rng.choice(np.array([-1., 1.]), size=(count, a["T"]))
                numerator = (xi*a["bar_h"]) @ a["e"]
                projected = np.zeros((count, a["N"]))
                for col in range(a["Q"].shape[1]):
                    projected += ((xi*a["Q"][:, col]) @ a["e"])**2
                variance = (ss-projected)/a["v"]
                if variance.min() <= 0:
                    raise ValueError("Nonpositive bootstrap residual variance")
                refitted[start:start+count] = core.ordered_sums(numerator**2/variance, order)
            for mode, values in (("frozen", frozen), ("refit", refitted)):
                result = core.evaluate(obs, values, a["N"], a["T"])
                rows.append(dict(setting=item["key"], rep=rep, denominator=mode, p_value=result["p"],
                                 reject_005=result["p"]<=.05, outer_seed=outer, inner_seed=inner))
                for j, k in enumerate(order):
                    comp.append(dict(setting=item["key"], rep=rep, denominator=mode, k=k,
                                     p_value=result["component_p"][j]))
            if (rep+1) % 50 == 0:
                print("Fitted-null %s: %d/%d" % (item["key"], rep+1, repetitions), flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(output_dir/"fitted_null_raw.csv", index=False)
    pd.DataFrame(comp).to_csv(output_dir/"fitted_null_components.csv", index=False)
    summary = []
    for (key, mode), g in frame.groupby(["setting", "denominator"]):
        rate, n = float(g.reject_005.mean()), len(g)
        z = 1.959963984540054
        center = (rate+z*z/(2*n))/(1+z*z/n)
        half = z*np.sqrt(rate*(1-rate)/n+z*z/(4*n*n))/(1+z*z/n)
        summary.append(dict(setting=key, denominator=mode, R=n, B_inner=inner_B, rejection_rate=rate,
                            rejections=int(g.reject_005.sum()), mcse=np.sqrt(rate*(1-rate)/n),
                            wilson_low=center-half, wilson_high=center+half))
    pd.DataFrame(summary).to_csv(output_dir/"fitted_null_summary.csv", index=False)


def break_reselection(item, output_dir, B):
    """Original trimmed scan, refitted bootstrap, and paired fixed-date comparison."""
    y, f = item["y"], item["base_f"]
    T, N = y.shape
    order = np.array(sorted({1, *core.grid(N)}), dtype=int)
    candidates = breaks_core.candidate_breaks(T, .2, 12)
    _, Q0, ops = breaks_core.score_operators(f, candidates)
    cache = breaks_core.fit_cache(f, candidates)
    observed_break, criterion = breaks_core.scan_break(y, Q0, ops, candidates)
    t2, coef, _ = breaks_core.fit_piecewise_t2(y, cache, observed_break)
    observed = core.ordered_sums(t2, order)
    fitted = cache[observed_break]["X"][:, 1:] @ coef[1:]
    residual = y-cache[observed_break]["X"] @ coef
    seed = 20260928+B+1000
    rng = np.random.default_rng(seed)
    selected = np.empty(B, dtype=int)
    boot, fixed = np.empty((B, len(order))), np.empty((B, len(order)))
    for start in range(0, B, 32):
        nb = min(32, B-start)
        xi = rng.choice(np.array([-1., 1.]), size=(nb, T))
        ystar = fitted[None, :, :]+xi[:, :, None]*residual[None, :, :]
        e0 = ystar-np.einsum("tr,rs,bsn->btn", Q0, Q0.T, ystar, optimize=True)
        s2 = np.maximum(np.sum(e0*e0, axis=1)/(T-Q0.shape[1]), 1e-20)
        score = np.einsum("cpt,btn->bcpn", ops, e0/np.sqrt(s2)[:, None, :], optimize=True)
        chosen = candidates[np.argmax(np.mean(np.sum(score*score, axis=2), axis=2), axis=1)]
        selected[start:start+nb] = chosen
        for j, t0 in enumerate(chosen):
            boot[start+j] = core.ordered_sums(breaks_core.fit_piecewise_t2(ystar[j], cache, int(t0))[0], order)
            fixed[start+j] = core.ordered_sums(breaks_core.fit_piecewise_t2(ystar[j], cache, observed_break)[0], order)
    rows, comp = [], []
    for name, values in (("reselected", boot), ("fixed_observed_break", fixed)):
        result = core.evaluate(observed[1:], values[:, 1:], N, T)
        qlo, qhi = np.quantile(selected, [.025, .975], method="nearest").astype(int)
        rows.append(dict(mode=name, B=B, seed=seed, observed_break_index=int(observed_break),
                         observed_break_month=str(item["months"][observed_break]),
                         observed_break_reselection_frequency=float(np.mean(selected==observed_break)),
                         selected_break_q025=str(item["months"][qlo]),
                         selected_break_q975=str(item["months"][qhi]),
                         AdaL_p=result["p"], cauchy_score=result["score"],
                         T1_p=(1+np.sum(values[:, 0]>=observed[0]))/(B+1)))
        probabilities = (1+np.sum(values>=observed, axis=0))/(B+1)
        for j, k in enumerate(order):
            comp.append(dict(mode=name, k=k, statistic=observed[j], p_value=probabilities[j]))
    pd.DataFrame(rows).to_csv(output_dir/"break_reselection_results.csv", index=False)
    pd.DataFrame(comp).to_csv(output_dir/"break_reselection_components.csv", index=False)
    values, count = np.unique(selected, return_counts=True)
    pd.DataFrame(dict(index=values, month=[str(item["months"][i]) for i in values], count=count,
                      probability=count/B)).to_csv(output_dir/"break_distribution.csv", index=False)
    pd.DataFrame(dict(index=candidates, month=[str(item["months"][i]) for i in candidates],
                      score=criterion)).to_csv(output_dir/"break_observed_scan.csv", index=False)


def plot_alphas(alphas, output_dir):
    cache = output_dir/"plot-cache"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(cache)
    os.environ["XDG_CACHE_HOME"] = str(cache)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["axes.unicode_minus"] = False
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(alphas.alpha_bp_monthly, bins=35, color="#2879a8", edgecolor="white")
    axes[0].set(xlabel="Monthly alpha (basis points)", ylabel="Number of assets")
    axes[1].hist(alphas.standardized_alpha, bins=35, color="#2879a8", edgecolor="white")
    for value in (-1.96, 1.96):
        axes[1].axvline(value, linestyle="--", color="#a43b36")
    axes[1].set(xlabel="Standardized alpha", ylabel="Number of assets")
    fig.tight_layout()
    fig.savefig(output_dir/"alpha_distribution.pdf")
    plt.close(fig)


def selected_paths(item, alphas, output_dir, plots):
    positions = np.argsort(-np.abs(alphas.standardized_alpha.to_numpy()), kind="stable")[:8]
    selected = alphas.iloc[positions].asset.astype(str).tolist()
    values = item["normalized"][selected]
    normalized = 100*values.div(values.iloc[0])
    normalized.to_csv(output_dir/"selected_stock_paths.csv", index_label="month")
    if plots:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 4))
        for asset in normalized:
            ax.plot(normalized.index.to_timestamp(), normalized[asset], label=asset, linewidth=1.1)
        ax.set(xlabel="Month", ylabel="Return/adjusted-price index (first month = 100)")
        ax.legend(ncol=4, fontsize=8)
        fig.tight_layout()
        fig.savefig(output_dir/"selected_stock_paths.pdf")
        plt.close(fig)


def input_manifest(data_dir):
    records = []
    paths = list(Path(data_dir).glob("China/CSMAR/*.csv")) + list(Path(data_dir).glob("America/*.csv")) + list(Path(data_dir).glob("America/*.CSV"))
    for path in sorted(set(paths)):
        h = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024*1024), b""):
                h.update(chunk)
        records.append(dict(file=str(path.relative_to(data_dir)), bytes=path.stat().st_size, sha256=h.hexdigest()))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("test", "reproduce", "check"))
    parser.add_argument("--data-dir")
    parser.add_argument("--output-dir")
    parser.add_argument("--y")
    parser.add_argument("--f")
    parser.add_argument("--analysis", default="all", choices=("all", "china", "us_ff3", "us_ff5", "break_reselection", "fitted_null", "diagnostics", "sensitivity"))
    parser.add_argument("--B", type=int)
    parser.add_argument("--B-maxk", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260929)
    parser.add_argument("--sensitivity-seed", type=int, default=20260906)
    parser.add_argument("--repetitions", type=int, default=500)
    parser.add_argument("--inner-B", type=int, default=999)
    parser.add_argument("--skip-maxk", action="store_true")
    parser.add_argument("--diagnostics", action="store_true")
    parser.add_argument("--plots", action="store_true")
    args = parser.parse_args()
    if args.command == "check":
        print(json.dumps(dict(python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__)))
        return
    for name in ("B", "B_maxk", "repetitions", "inner_B"):
        value = getattr(args, name)
        if value is not None and value < 1:
            parser.error(name+" must be positive")
    if args.seed < 0 or args.sensitivity_seed < 0:
        parser.error("seed and sensitivity-seed must be nonnegative")
    if not args.output_dir:
        parser.error("--output-dir is required")
    out = Path(args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    versions = dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
                    platform=platform.platform(), arguments=vars(args))
    with threadpool_limits(limits=1):
        if args.command == "test":
            if not args.y or not args.f:
                parser.error("test needs --y and --f CSV matrices without headers")
            y, f = np.loadtxt(args.y, delimiter=",", ndmin=2), np.loadtxt(args.f, delimiter=",", ndmin=2)
            empirical_test(y, f, args.B or 9999, args.seed, not args.skip_maxk, args.B_maxk, out)
        else:
            if not args.data_dir:
                parser.error("reproduce needs --data-dir")
            requested = {"all": ["china", "us_ff3", "us_ff5"], "fitted_null": ["china", "us_ff3"],
                         "diagnostics": ["china", "us_ff3"], "sensitivity": ["china", "us_ff3"],
                         "break_reselection": ["us_ff3"]}.get(args.analysis, [args.analysis])
            items = load_settings(args.data_dir, out, requested)
            versions["inputs"] = input_manifest(args.data_dir)
            if args.analysis not in ("fitted_null", "break_reselection", "diagnostics", "sensitivity"):
                rows, components = [], []
                for item in items:
                    dest = out/item["key"]
                    table, component, alphas = empirical_test(item["y"], item["f"], args.B or item["B"],
                        item["seed"], not args.skip_maxk, args.B_maxk,
                        dest, item["assets"], item["key"], item["maxk_seed"], item["consume_l_before_maxk"])
                    rows.append(table)
                    component.insert(0, "setting", item["key"])
                    components.append(component)
                    if args.diagnostics:
                        residual_diagnostics(item, dest)
                        model_diagnostics(item, dest, args.B or (4999 if item["key"]=="us_ff5" else 1999))
                    if args.plots:
                        plot_alphas(alphas, dest)
                    selected_paths(item, alphas, dest, args.plots)
                    magnitude = alphas.alpha_bp_monthly
                    pd.DataFrame([dict(setting=item["key"], median_bp=magnitude.median(),
                        q05_bp=magnitude.quantile(.05), q95_bp=magnitude.quantile(.95),
                        mean_abs_bp=magnitude.abs().mean(), max_abs_bp=magnitude.abs().max(),
                        share_abs_ge_50bp=(magnitude.abs()>=50).mean(),
                        BH_005=int((alphas.BH_q<=.05).sum()), BH_010=int((alphas.BH_q<=.10).sum()))]
                        ).to_csv(dest/"economic_magnitudes.csv", index=False)
                pd.concat(rows, ignore_index=True).to_csv(out/"main_results.csv", index=False)
                pd.concat(components, ignore_index=True).to_csv(out/"component_probabilities.csv", index=False)
            if args.analysis == "diagnostics":
                for item in items:
                    dest = out/item["key"]
                    dest.mkdir(parents=True, exist_ok=True)
                    residual_diagnostics(item, dest)
                    model_diagnostics(item, dest, args.B or 1999)
            if args.analysis in ("all", "sensitivity"):
                if args.skip_maxk and args.analysis == "sensitivity":
                    raise ValueError("analysis='sensitivity' requires include_maxk=TRUE")
                if not args.skip_maxk:
                    penalty_sensitivity(items, out, args.B_maxk, args.sensitivity_seed)
                else:
                    versions["maxk_sensitivity"] = "omitted because include_maxk=FALSE"
            if args.analysis in ("all", "fitted_null"):
                fitted_null(items, out, args.repetitions, args.inner_B)
            if args.analysis in ("all", "break_reselection"):
                break_reselection(next(x for x in items if x["key"]=="us_ff3"), out, args.B or 9999)
    save_json(versions, out/"provenance.json")
    print("Completed empirical workflow: "+str(out), flush=True)


if __name__ == "__main__":
    main()
