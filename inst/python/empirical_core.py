"""Numerical functions extracted from the dated research scripts.
See SOURCE_MAP.json for provenance; no market observations are distributed.
"""
from __future__ import annotations
import hashlib
import numpy as np
from scipy import stats
from quic_solver import solve as quic_solve
BASESEED = 20260903
def grid(n):
    if n < 40:
        raise ValueError("Ma grid requires N>=40")
    count = int(np.floor(np.log2(n / 20)))
    result = np.array(sorted({5, *(int(np.ceil(n / 2**j)) for j in range(1, count + 1))}), dtype=int)
    assert len(result) == count + 1
    return result

def fit(y, f):
    """Use the original saved analysis's OLS arithmetic for reproducibility."""
    x = np.column_stack([np.ones(len(f)), np.asarray(f)])
    qx, _ = np.linalg.qr(x)
    inv = np.linalg.inv(x.T @ x)
    coef = inv @ x.T @ np.asarray(y)
    residual = np.asarray(y) - x @ coef
    df = len(x) - x.shape[1]
    variance = np.sum(residual * residual, axis=0) / df
    info = 1 / inv[0, 0]
    w = inv[0] @ x.T
    qf, _ = np.linalg.qr(x[:, 1:], mode="reduced")
    mf_one = np.ones(len(x)) - qf @ (qf.T @ np.ones(len(x)))
    effective = mf_one @ mf_one
    bar_h = mf_one / np.sqrt(effective)
    np.testing.assert_allclose(bar_h, np.sqrt(info) * w, rtol=0, atol=5e-12)
    np.testing.assert_allclose(effective, info, rtol=0, atol=5e-10)
    assert np.min(variance) > 0 and df > 0
    return dict(X=x, Q=qx, coef=coef, e=residual, v=df, s2=variance, q=info, w=w,
                bar_h=bar_h, t=coef[0] * np.sqrt(info / variance), T=len(x), N=residual.shape[1])

def ordered_sums(t2, orders):
    return np.cumsum(np.sort(t2, axis=-1)[..., ::-1], axis=-1)[..., orders - 1]

def bootstrap(a, orders, draws, seed, batch=256):
    rng = np.random.default_rng(seed)
    result = np.empty((draws, len(orders)))
    for start in range(0, draws, batch):
        n = min(batch, draws - start)
        xi = rng.choice(np.array([-1.0, 1.0]), size=(n, a["T"]))
        numerator = (xi * a["bar_h"]) @ a["e"]
        t2 = numerator**2 / a["s2"]
        result[start:start + n] = ordered_sums(t2, orders)
    return result

def evaluate(observed, draws, n, t):
    exceedance = np.sum(draws >= observed[None, :], axis=0)
    component = (1 + exceedance) / (len(draws) + 1)
    epsilon = 1 / (n + t)
    clipped = np.clip(component, epsilon, 1 - epsilon)
    score = np.mean(np.tan(np.pi * (0.5 - clipped)))
    probability = 0.5 - np.arctan(score) / np.pi
    return dict(exceedance=exceedance, component_p=component, clipped_p=clipped,
                epsilon=epsilon, score=float(score), p=float(probability), reject_005=bool(probability <= 0.05))
def seed_from_legacy(base: int, key: str) -> int:
    """Exact deterministic seed rule used by the archived empirical scripts."""
    return (base + int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)) % (2**32)

def original_stream_bootstrap(a: dict, orders: np.ndarray, draws: int, seed: int):
    """Extend the old refitted-denominator bootstrap to additional top-k sums.

    The arithmetic and batch size intentionally match
    ``real_data_20260903.l_bootstrap`` so archived common columns can be checked
    draw by draw.  The legacy 1e-20 safeguard is retained but its activation is
    audited; it is not expected to bind in these samples.
    """
    rng = np.random.default_rng(seed)
    result = np.empty((draws, len(orders)))
    residual = a["e"]
    total_ss = np.sum(residual * residual, axis=0)
    minimum_variance = float("inf")
    safeguard_count = 0
    batch_size = 128
    for start in range(0, draws, batch_size):
        count = min(batch_size, draws - start)
        xi = rng.choice(np.array([-1.0, 1.0]), size=(count, a["T"]))
        alpha_star = (xi * a["w"]) @ residual
        projected_ss = np.zeros((count, a["N"]))
        for column in range(a["Q"].shape[1]):
            projected_ss += ((xi * a["Q"][:, column]) @ residual) ** 2
        variance = (total_ss - projected_ss) / a["v"]
        minimum_variance = min(minimum_variance, float(np.min(variance)))
        safeguard_count += int(np.sum(variance <= 1e-20))
        variance = np.maximum(variance, 1e-20)
        t_squared = a["q"] * alpha_star**2 / variance
        result[start : start + count] = ordered_sums(t_squared, orders)
    return result, {
        "batch_size": batch_size,
        "minimum_unclipped_refitted_variance": minimum_variance,
        "legacy_variance_safeguard_count": safeguard_count,
    }
def seed(key): return (BASESEED+int(hashlib.sha256(key.encode()).hexdigest()[:8],16))%(2**32)

def py_test(a,t2=None):
    t2=a['t']**2 if t2 is None else t2;N=a['N'];v=a['v']
    e=a['e'];z=e/np.sqrt(np.sum(e*e,axis=0));corr=z.T@z
    off=corr[np.triu_indices(N,1)];th=stats.norm.isf(.1/(2*(N-1)))**2
    rho2=np.mean(off**2*(v*off**2>=th));cent=v/(v-2)
    denom=cent*np.sqrt(2*(v-1)/(v-4)*(1+(N-1)*rho2))
    val=(np.sum(t2)-N*cent)/np.sqrt(N)/denom
    return float(val),float(stats.norm.sf(val)),float(rho2),float(denom)

def cauchy(p):
    # Finite bootstrap ranks use the midpoint grid, excluding exactly 0 and 1.
    p=np.asarray(p);s=np.mean(1/np.tan(np.pi*p),axis=-1)
    return s,stats.cauchy.sf(s)

def l_values(t2):
    N=t2.shape[-1];k=[1,5,int(np.ceil(.25*N)),int(np.ceil(.5*N))]
    cs=np.cumsum(np.sort(t2,axis=-1)[...,::-1],axis=-1)
    return cs[...,np.array(k)-1]

def l_bootstrap(a,B,rng,block=1):
    e=a['e'];T,N=e.shape;out=[];Q=a['Q'];tot=np.sum(e*e,axis=0)
    for start in range(0,B,128):
        n=min(128,B-start)
        if block==1: xi=rng.choice([-1.,1.],size=(n,T))
        else:
            # Randomly shifted nonoverlapping blocks; preserves within-block temporal products.
            xi=np.empty((n,T))
            for b in range(n):
                offset=rng.integers(block);labs=(np.arange(T)+offset)//block
                xi[b]=rng.choice([-1.,1.],size=labs.max()+1)[labs]
        ast=(xi*a['w'])@e
        proj=np.zeros((n,N))
        for j in range(Q.shape[1]):proj+=((xi*Q[:,j])@e)**2
        s2=np.maximum((tot-proj)/a['v'],1e-20)
        out.append(l_values(a['q']*ast**2/s2))
    return np.concatenate(out)

def precision(a,rho_scale=1.):
    T,N=a['T'],a['N'];delta=np.log(np.log(T))*np.sqrt(np.log(N))
    se=np.sqrt(a['s2']*a['v']/T/a['q'])
    screen=np.abs(a['coef'][0])/se>delta
    U=a['e']+a['coef'][0][None,:]*(~screen)[None,:]
    cov=U.T@U/T;sd=np.sqrt(np.diag(cov));R=cov/sd[:,None]/sd[None,:]
    rho=rho_scale*np.sqrt(np.log(N)/T)
    cv,prec,solver_meta=quic_solve(R,rho,warm=True)
    G=prec/sd[:,None]/sd[None,:]
    return G,dict(glasso_rho=rho,**solver_meta,
                   precision_offdiag_density=float(np.mean(np.abs(prec[np.triu_indices(N,1)])>1e-10)),screen_count=int(screen.sum()))

def xia_test(a,B,rng,rho_scale=1.):
    if B < 2:
        raise ValueError('max-k requires at least two bootstrap draws')
    G,meta=precision(a,rho_scale);T,N=a['T'],a['N'];K=min(30,N)
    obs=T*(G@a['coef'][0])**2/np.diag(G);go=np.cumsum(np.sort(obs)[::-1])[:K]
    trans=a['e']@G/np.sqrt(np.diag(G))[None,:]/np.sqrt(a['q'])
    boot=[]
    for j in range(0,B,128):
        v=rng.standard_normal((min(128,B-j),T))@trans
        boot.append(np.cumsum(np.sort(v*v,axis=1)[:,::-1],axis=1)[:,:K])
    gb=np.concatenate(boot);mu=gb.mean(0);sd=gb.std(0,ddof=0)
    if not np.isfinite(sd).all() or np.any(sd <= 0):
        raise ValueError('max-k bootstrap standard deviations must be finite and positive')
    zo=(go-mu)/sd;zb=(gb-mu)/sd
    rows={}
    for maxk in [5,10,30]:
        k=min(maxk,K);stat=float(zo[:k].max());pv=(1+np.sum(zb[:,:k].max(1)>=stat))/(B+1)
        rows[f'maxk_K{maxk}_stat']=stat;rows[f'maxk_K{maxk}_p']=float(pv);rows[f'maxk_K{maxk}_selected_k']=int(np.argmax(zo[:k])+1)
    return rows,meta,G,gb
