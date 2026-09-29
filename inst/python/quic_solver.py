"""Lazy portable QUIC binding; see vendor/quic/SOURCE.md and LICENSE."""
import ctypes as c
from pathlib import Path
import numpy as np

_fun=None
_cache_dir=None

def configure(cache_dir):
    global _cache_dir
    _cache_dir=Path(cache_dir)

def _load():
    global _fun
    if _fun is not None:
        return
    import os, platform, shutil, subprocess, sys
    if _cache_dir is None:
        raise RuntimeError("QUIC cache directory is not configured")
    system=platform.system()
    if system not in ("Darwin","Linux"):
        raise RuntimeError("The empirical QUIC backend currently supports macOS and Linux; use include_maxk=FALSE on Windows")
    compiler=shutil.which(os.environ.get("CXX","c++"))
    if not compiler:
        raise RuntimeError("max-k requires a C++ compiler and BLAS/LAPACK; use include_maxk=FALSE to omit it")
    _cache_dir.mkdir(parents=True,exist_ok=True)
    library=_cache_dir/("libquic.dylib" if system=="Darwin" else "libquic.so")
    source=Path(__file__).parent/'vendor/quic/QUIC.C'
    if not library.exists():
        flags=(["-dynamiclib","-framework","Accelerate"] if system=="Darwin" else ["-shared","-fPIC","-llapack","-lblas"])
        command=[compiler,"-O3",str(source),"-o",str(library)]+flags
        done=subprocess.run(command,capture_output=True,text=True)
        if done.returncode:
            raise RuntimeError("QUIC compilation failed: "+done.stderr)
    lib=c.CDLL(str(library))
    double=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    uint=np.ctypeslib.ndpointer(dtype=np.uint32,flags='C_CONTIGUOUS')
    _fun=lib.QUIC
    _fun.argtypes=[c.c_char,c.POINTER(c.c_uint32),double,double,c.POINTER(c.c_uint32),double,
                  c.POINTER(c.c_double),c.POINTER(c.c_int32),c.POINTER(c.c_uint32),
                  double,double,double,double,uint,double]
    _fun.restype=None

_warm={}

def solve(R,rho,tol=1e-7,max_iter=200,warm=False):
    _load()
    R=np.ascontiguousarray(R,dtype=float);n=R.shape[0]
    lam=np.full((n,n),rho);np.fill_diagonal(lam,0)
    key=(n,float(rho))
    if warm and key in _warm:
        P,W=[v.copy() for v in _warm[key]]
    else:P=np.eye(n);W=np.eye(n)
    opt=np.zeros(1);cpu=np.zeros(1);gap=np.zeros(1);it=np.zeros(1,dtype=np.uint32)
    nn=c.c_uint32(n);pl=c.c_uint32(1);tt=c.c_double(tol);msg=c.c_int32(0);mi=c.c_uint32(max_iter)
    _fun(b'd',c.byref(nn),R,lam,c.byref(pl),np.ones(1),c.byref(tt),c.byref(msg),c.byref(mi),P,W,opt,cpu,it,gap)
    # KKT residual uses W=P^{-1}, which QUIC maintains explicitly.
    nz=np.abs(P)>1e-9;off=~np.eye(n,dtype=bool)
    err=np.zeros_like(R);err[nz&off]=np.abs((R-W+lam*np.sign(P))[nz&off])
    err[~nz&off]=np.maximum(np.abs((R-W)[~nz&off])-rho,0)
    err[np.diag_indices(n)]=np.abs(np.diag(R-W))
    kkt=float(err.max())
    if it[0]>=max_iter or not np.isfinite(P).all() or kkt>1e-4:
        raise RuntimeError(f'QUIC did not converge: iterations={it[0]}, gap={gap[0]}, KKT={kkt}')
    if warm:
        _warm.clear();_warm[key]=(P.copy(),W.copy())
    return W,P,dict(solver='QUIC 1.2',glasso_iterations=int(it[0]),glasso_dual_gap=float(gap[0]),glasso_kkt_max=kkt,glasso_objective=float(opt[0]))
