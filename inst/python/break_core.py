"""Pure score-scan and OLS routines from the final break-reselection analysis."""
from __future__ import annotations
import numpy as np
def candidate_breaks(T: int, trim: float, min_segment: int) -> np.ndarray:
    """Post-break start indices, matching core.pretests inclusive endpoints."""
    lo = max(min_segment, int(trim * T))
    hi = min(T - min_segment, int((1.0 - trim) * T))
    if lo > hi:
        raise ValueError("trim/min_segment leave no admissible break date")
    return np.arange(lo, hi + 1, dtype=int)

def score_operators(f: np.ndarray, breaks: np.ndarray):
    """Operators for the common slope-change score scan."""
    T = len(f)
    X0 = np.column_stack([np.ones(T), f])
    Q0, _ = np.linalg.qr(X0)
    ops = []
    for t0 in breaks:
        h = f * (np.arange(T) >= t0)[:, None]
        h = h - Q0 @ (Q0.T @ h)
        qh, rh = np.linalg.qr(h)
        if np.min(np.abs(np.diag(rh))) < 1e-10:
            raise np.linalg.LinAlgError(f"rank-deficient break regressors at {t0}")
        ops.append(qh.T)
    return X0, Q0, np.asarray(ops)

def scan_break(y: np.ndarray, Q0: np.ndarray, ops: np.ndarray, breaks: np.ndarray):
    """Select the break maximizing the cross-sectional mean score norm."""
    T = y.shape[0]
    e0 = y - Q0 @ (Q0.T @ y)
    s2 = np.maximum(np.sum(e0 * e0, axis=0) / (T - Q0.shape[1]), 1e-20)
    es = e0 / np.sqrt(s2)
    score = np.einsum("cpt,tn->cpn", ops, es, optimize=True)
    criterion = np.mean(np.sum(score * score, axis=1), axis=1)
    pos = int(np.argmax(criterion))
    return int(breaks[pos]), criterion

def piecewise_design(f: np.ndarray, t0: int, intercept: bool = True) -> np.ndarray:
    post = (np.arange(len(f)) >= t0).astype(float)
    cols = [f, f * post[:, None]]
    if intercept:
        cols.insert(0, np.ones((len(f), 1)))
    return np.column_stack(cols)

def fit_cache(f: np.ndarray, breaks: np.ndarray):
    """Cache OLS maps for all candidate piecewise models."""
    out = {}
    for t0 in breaks:
        X = piecewise_design(f, int(t0), intercept=True)
        inv = np.linalg.inv(X.T @ X)
        out[int(t0)] = {
            "X": X,
            "inv": inv,
            "w": inv[0] @ X.T,
            "q": 1.0 / inv[0, 0],
            "v": len(f) - X.shape[1],
        }
    return out

def fit_piecewise_t2(y: np.ndarray, cache: dict, t0: int):
    c = cache[int(t0)]
    xty = c["X"].T @ y
    coef = c["inv"] @ xty
    sse = np.maximum(np.sum(y * y, axis=0) - np.sum(coef * xty, axis=0), 1e-20)
    s2 = sse / c["v"]
    t2 = c["q"] * coef[0] ** 2 / s2
    return t2, coef, s2
