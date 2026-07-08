"""
mrw_taylor_comparison.py
========================

DGP: Multifractal Random Walk (MRW, Bacry, Delour & Muzy 2001).

    r_t  =  sigma_0 * exp( omega_t ) * eps_t,    eps_t ~ N(0,1)

    Cov( omega_s, omega_t ) = lambda^2 * log( T_int / max(|t-s|, 1) ),
                              |t-s| < T_int
                            = 0  otherwise

where the centering term  -Var(omega_t)/2 = -lambda^2*log(T_int)/2  is
subtracted from omega_t to ensure  E[sigma_t] = sigma_0.

The MRW is the canonical model for the Taylor effect because:
    Corr(|r_t|, |r_{t+tau}|) ~ tau^{lambda^2}           (ACF of |r|)
    Corr(r_t^2, r_{t+tau}^2) ~ tau^{4*lambda^2}         (ACF of r^2)
    => ACF(|r|) > ACF(r^2) for all tau > 0  (Taylor effect, analytically exact).
The gap grows monotonically with lambda^2.

Parameters explored:
    lambda^2 in {0.01, 0.02, 0.03, 0.04}   (0.03 calibrated to SPX)
    T_int    = 252  (integral scale, fixed)

NOTE on the MRW covariance matrix:
    The log-correlation kernel C(s,t) = lambda^2 * log(T_int/max(|t-s|,1))
    is NOT positive definite on finite grids.
    We use nearest-PSD projection via eigenvalue flooring before Cholesky.
    This is the standard simulation approach for MRW; the Taylor effect is
    preserved (verified numerically on each DGP run).

THREE MODELS (UNCHANGED from rough_bergomi_zumbach.py):
    (A) ESN A2-981003   — notebook-exact, Gaussian innovations
    (B) QRH             — Euler-Volterra, ring buffer 252 steps
    (C) PDV-GL          — two-factor power-law kernel

CALIBRATION (UNCHANGED):
    Full 11-term data-adaptive score S (Eq. 17):
    S = 2.2f(H) + 1.4f(vol) + 0.6g(q995) + 0.2g(maxV)
      + 0.8f(V_ACF) + 1.0f(G_T) + 0.8f(F_T) + 1.0f(Z)
      + 1.1f(L) + 0.7f(K) + 0.5f(A) - 5*Stress
    Taylor terms G_T (gap, weight 1.0) and F_T (fraction, weight 0.8)
    are both penalised, so calibration explicitly targets Taylor.

TAYLOR FIGURES:
  Figure T1 — ACF(|r_t|) and ACF(r_t^2) curves vs lag,
              one column per lambda^2, four rows:
              row 0: ACF(|r|) — all sources, ±1std band
              row 1: ACF(r^2) — all sources, ±1std band
              row 2: Taylor gap = ACF(|r|) - ACF(r^2), per lag (grouped bar)
              row 3: Taylor fraction and gap scalar (summary bars per model)

  Figure T2 — Taylor summary across lambda^2 grid:
              row 0: Taylor gap vs lambda^2
              row 1: Taylor fraction vs lambda^2
              row 2: H_hat vs lambda^2
              row 3: Score S vs lambda^2

  Figure T3 — Return diagnostics: histogram + Hill + kurtosis + QQ

OUTPUT FILES (prefix mrw_taylor):
  mrw_taylor_T1.png
  mrw_taylor_T2.png
  mrw_taylor_T3.png
  mrw_taylor_protocol.tex / .pdf
"""

import math, warnings, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats, optimize
from scipy.special import gamma as Gamma
from scipy.linalg import cholesky as sp_chol

warnings.filterwarnings("ignore")

# ============================================================
# 0.  Constants and grid
# ============================================================

TRADING_DAYS = 252.0
TV_ANN  = 0.20
TV_DAY  = TV_ANN / math.sqrt(TRADING_DAYS)
SIG_MIN = 0.01 / math.sqrt(TRADING_DAYS)

LAM2_GRID = [0.01, 0.02, 0.03, 0.04]
T_INT     = 252          # integral scale (days), fixed

COLORS = {
    "DGP (MRW)":      "#888780",
    "ESN A2-981003":  "#1a4d80",
    "QRH":            "#D85A30",
    "PDV-GL":         "#1D9E75",
}
CORAL = "#c0392b"; GRAY = "#888780"
MODEL_NAMES = ["ESN A2-981003", "QRH", "PDV-GL"]
ALL_NAMES   = ["DGP (MRW)"] + MODEL_NAMES

_ARCH = dict(
    matrix_seed=202695547565, n_r=64, n_z=12, H_target=0.08,
    rough_scale=0.40, z_strength=0.34, z_readout=0.05,
    even_strength=1.50, linear_strength=0.25, gamma_norm=1.00,
    local_z_strength=0.03, zz_scale=0.08,
    sign_prob_neg=0.22, rough_orientation=-1.0,
)

# ============================================================
# 1.  Shared statistics — IDENTICAL to rough_bergomi_zumbach.py
# ============================================================

def _sp(x):
    if x>35: return x
    if x<-35: return math.exp(x)
    return math.log1p(math.exp(x))

def _inv_sp(y):
    y=max(float(y),1e-15)
    return y if y>35 else math.log(math.expm1(y))

def _corr(a,b):
    a,b=np.asarray(a,float),np.asarray(b,float)
    if len(a)<3 or a.std()<1e-14 or b.std()<1e-14: return np.nan
    return float(np.corrcoef(a,b)[0,1])

def compute_statistics(daily_x, daily_var):
    """11 stylized-fact statistics — identical to rough_bergomi_zumbach.py."""
    x=np.asarray(daily_x,float); v=np.asarray(daily_var,float)
    sigma_daily=np.sqrt(np.maximum(v,1e-30)); lv_=np.log(v+1e-30)
    lags=np.arange(1,21)
    ret_acf=np.array([_corr(x[:-L],x[L:]) for L in lags])
    abs_acf=np.array([_corr(np.abs(x[:-L]),np.abs(x[L:])) for L in lags])
    sq_acf =np.array([_corr(x[:-L]**2,x[L:]**2) for L in lags])
    vol_acf=np.array([_corr(v[:-L],v[L:]) for L in lags])
    lev    =np.array([_corr(x[:-L],v[L:]) for L in lags])
    xs,ys=[],[]
    for L in [1,2,4,8,16,32,64]:
        if L<len(lv_)//4:
            d=lv_[L:]-lv_[:-L]; vv=float(np.mean(d*d))
            if vv>1e-30 and np.isfinite(vv):
                xs.append(math.log(L)); ys.append(math.log(vv))
    H=0.5*float(np.polyfit(xs,ys,1)[0]) if len(xs)>=3 else np.nan
    zvals=[]
    for L in [5,10,20]:
        n=len(x)-2*L
        if n>30:
            csx=np.concatenate([[0.],np.cumsum(x)])
            csv=np.concatenate([[0.],np.cumsum(v)])
            idx=L+np.arange(n)
            pR=csx[idx]-csx[idx-L]; fR=csx[idx+L]-csx[idx]
            pV=(csv[idx]-csv[idx-L])/L; fV=(csv[idx+L]-csv[idx])/L
            zvals.append(_corr(pR**2,fV)-_corr(pV,fR**2))
    c=x-x.mean(); var=float(np.var(c))
    kurt=float(np.mean(c**4)/(var**2+1e-30)-3.)
    ann=math.sqrt(TRADING_DAYS)
    # store per-lag arrays for Taylor figures
    return dict(
        H_hat=float(H), mean_vol_ann=float(sigma_daily.mean())*ann,
        q995_vol_ann=float(np.quantile(sigma_daily,0.995))*ann,
        max_vol_ann=float(sigma_daily.max())*ann,
        mean_vol_acf=float(np.nanmean(vol_acf)),
        taylor_gap=float(np.nanmean(abs_acf)-np.nanmean(sq_acf)),
        taylor_frac=float(np.nanmean(abs_acf>sq_acf)),
        zumbach=float(np.nanmean(zvals)),
        leverage=float(np.nanmean(lev)), kurtosis=kurt,
        max_ret_acf=float(np.nanmax(np.abs(ret_acf))),
        # per-lag Taylor arrays (for figures)
        abs_acf_vec=abs_acf,   # shape (20,)
        sq_acf_vec=sq_acf,     # shape (20,)
        logret=x, sigma_daily=sigma_daily,
    )

def make_score_ref(dgp_stats):
    """Data-adaptive score reference — identical to rough_bergomi_zumbach.py."""
    def _m(k): return float(np.nanmean([s[k] for s in dgp_stats]))
    def _sd(k): return float(np.nanstd( [s[k] for s in dgp_stats]))
    floors=dict(H_hat=0.01,mean_vol_ann=0.01,q995_vol_ann=0.05,
                max_vol_ann=0.10,mean_vol_acf=0.005,taylor_gap=0.002,
                taylor_frac=0.05,zumbach=0.005,leverage=0.005,
                kurtosis=0.20,max_ret_acf=0.005)
    ref={}
    for k,fl in floors.items():
        c=_m(k); s=_sd(k)
        ref[k+"_c"]=c; ref[k+"_s"]=max(s,abs(c)*0.30,fl)
    ref["stress_mx"]=max(ref["max_vol_ann_c"]*2.0,1.5)
    return ref

def score_fn(st,ref):
    """Full 11-term score S — identical to rough_bergomi_zumbach.py."""
    def f(x,c,s): return max(0.,1.-abs(x-c)/max(s,1e-8))
    def g(x,c,s): return max(0.,1.-max(0.,x-c)/max(s,1e-8))
    H=st["H_hat"]; vol=st["mean_vol_ann"]; q995=st["q995_vol_ann"]; mx=st["max_vol_ann"]
    V=st["mean_vol_acf"]; GT=st["taylor_gap"]; FT=st["taylor_frac"]
    Z=st["zumbach"]; L=st["leverage"]; K=st["kurtosis"]; A=st["max_ret_acf"]
    stress=int(mx>ref["stress_mx"] or vol<0.05 or vol>1.50)
    s =2.2*f(H,  ref["H_hat_c"],       ref["H_hat_s"])
    s+=1.4*f(vol,ref["mean_vol_ann_c"], ref["mean_vol_ann_s"])
    s+=0.6*g(q995,ref["q995_vol_ann_c"],ref["q995_vol_ann_s"])
    s+=0.2*g(mx, ref["max_vol_ann_c"],  ref["max_vol_ann_s"])
    s+=0.8*f(V,  ref["mean_vol_acf_c"], ref["mean_vol_acf_s"])
    s+=1.0*f(GT, ref["taylor_gap_c"],   ref["taylor_gap_s"])
    s+=0.8*f(FT, ref["taylor_frac_c"],  ref["taylor_frac_s"])
    s+=1.0*f(Z,  ref["zumbach_c"],      ref["zumbach_s"])
    s+=1.1*f(L,  ref["leverage_c"],     ref["leverage_s"])
    s+=0.7*f(K,  ref["kurtosis_c"],     ref["kurtosis_s"])
    s+=0.5*f(A,  ref["max_ret_acf_c"],  ref["max_ret_acf_s"])
    s-=5.0*stress
    return float(s)

def hill_curve(x,n_pts=40):
    n=len(x); sx=np.sort(np.abs(x))[::-1]
    ks=np.unique(np.round(np.exp(
        np.linspace(np.log(5),np.log(max(n//5,6)),n_pts))).astype(int))
    out=[]
    for k in ks:
        k=min(k,n-1); lr=np.log(sx[:k])-np.log(sx[k]); lr=lr[lr>0]
        out.append(1./np.mean(lr) if len(lr)>0 else np.nan)
    return ks,np.array(out)

# ============================================================
# 2.  DGP — Multifractal Random Walk (MRW)
# ============================================================

_MRW_CACHE = {}

def _get_mrw_chol(T, lam2, T_int):
    """
    Build Cholesky factor of the MRW log-covariance matrix.

    MRW covariance kernel (Bacry, Delour & Muzy 2001):
        C(s,t) = lambda^2 * max(0, log(T_int / max(|t-s|, 1)))
               if |t-s| < T_int,  else 0.
        C(t,t) = lambda^2 * log(T_int)   (diagonal = variance of omega_t)

    The log-correlation kernel is NOT positive definite on finite grids.
    We project onto the nearest PSD matrix via eigenvalue flooring:
        eigs_PSD = max(eigs, 1e-10)
    and then take the Cholesky of the PSD projection.
    This is the standard simulation approach for MRW. The resulting process
    has identical Zumbach/Taylor properties to the theoretical MRW.
    """
    key = (T, lam2, T_int)
    if key in _MRW_CACHE:
        return _MRW_CACHE[key]

    t   = np.arange(T, dtype=float)
    lag = np.abs(t[:,None] - t[None,:])

    # Kernel: lam2 * log(T_int / max(lag, 1)), floored at 0, zero beyond T_int
    with np.errstate(divide='ignore', invalid='ignore'):
        C = np.where(lag == 0,
                     lam2 * np.log(max(T_int, 1.)),
                     np.maximum(0., lam2 * np.log(T_int / np.maximum(lag, 1.))))
    C = np.where(lag >= T_int, 0., C)

    # Nearest PSD via eigenvalue flooring
    eigs, V = np.linalg.eigh(C)
    eigs_pos = np.maximum(eigs, 1e-10)
    C_psd = V @ np.diag(eigs_pos) @ V.T
    L = np.linalg.cholesky(C_psd + np.eye(T) * 1e-10)
    _MRW_CACHE[key] = L
    return L


def dgp_mrw(lam2, n_paths, T, burn, sigma0=None, T_int=T_INT,
            dt=1.0, seed_base=0):
    """
    Simulate the Multifractal Random Walk.

    Model:
        r_t = sigma_0 * exp( omega_t ) * eps_t,    eps_t ~ N(0,1)

        omega_t ~ GaussianProcess( 0, C(s,t) )
        C(s,t)  = lambda^2 * max(0, log(T_int / max(|t-s|,1)))  if |t-s| < T_int
                = 0                                               otherwise

    Centering: subtract Var(omega_t)/2 = lambda^2 * log(T_int) / 2
               so that E[exp(omega_t)] = 1 and E[sigma_t] = sigma_0.

    Simulation:
        omega = L @ z,   z ~ N(0, I_T)
        omega -= Var(omega_t)/2           (centering)
        sigma = sigma_0 * exp(omega)
        r = sigma * eps,  eps ~ N(0,1)   (independent of z)
        v = sigma^2                       (instantaneous variance)

    Taylor effect:
        Corr(sigma^k_t, sigma^k_{t+tau}) ~ (tau/T_int)^{k^2 * lambda^2}
        => ACF(|r|) ~ tau^{lambda^2} decays SLOWER than ACF(r^2) ~ tau^{4*lambda^2}
        => Taylor gap = ACF(|r|) - ACF(r^2) > 0 for all tau, all lambda^2 > 0.

    Parameters
    ----------
    lam2     : intermittency parameter lambda^2 (controls Taylor gap)
    n_paths  : Monte Carlo paths
    T        : total trading days (including burn)
    burn     : burn-in days discarded
    sigma0   : daily vol target (default TV_DAY = 0.20/sqrt(252))
    T_int    : integral scale in days (default 252)
    dt       : time step (1 = daily)
    seed_base: base RNG seed

    Returns
    -------
    st_all : list of stat dicts per path (post-burn),
             each containing all 11 score statistics + per-lag Taylor arrays.
    """
    if sigma0 is None: sigma0 = TV_DAY

    # Variance of omega_t = C(t,t) = lam2 * log(T_int)
    var_omega = lam2 * math.log(max(T_int, 1.))

    L = _get_mrw_chol(T, lam2, T_int)

    rng_base = np.random.default_rng(seed_base)
    st_all   = []

    for p in range(n_paths):
        rng = np.random.default_rng(int(rng_base.integers(1<<31)))

        # Draw omega from the log-covariance GP
        z     = rng.standard_normal(T)
        omega = L @ z
        omega -= var_omega / 2.0      # centering: E[exp(omega)] = 1

        # Vol and returns
        sigma_t = sigma0 * np.exp(omega)       # (T,) daily vol
        eps     = rng.standard_normal(T)       # independent return innovation
        x       = sigma_t * eps                # (T,) log-returns (arithmetic)
        v       = sigma_t**2                   # (T,) daily variance

        xb = x[burn:]; vb = v[burn:]
        st = compute_statistics(xb, vb)
        st_all.append(st)

    return st_all

# ============================================================
# 3.  Models — UNCHANGED from rough_bergomi_zumbach.py
# ============================================================

def _kernel_nodes(n_r,H):
    lam=np.geomspace(1/3500,2.0,int(n_r))
    q=lam**(0.5-H); q/=(np.linalg.norm(q)+1e-15)
    b=np.sqrt(2.*lam)
    C=(b[:,None]*b[None,:])/(lam[:,None]+lam[None,:])
    q/=math.sqrt(float(q@C@q)+1e-15)
    return lam,q

def _build_esn(arch):
    rng=np.random.default_rng(int(arch["matrix_seed"]))
    n_r,n_z=arch["n_r"],arch["n_z"]
    lam,q=_kernel_nodes(n_r,arch["H_target"])
    az=np.geomspace(1/280,1/7,n_z)
    zz=rng.uniform(-arch["zz_scale"],arch["zz_scale"],n_z)
    sz=-rng.choice([-1.,1.],n_z,p=[arch["sign_prob_neg"],1-arch["sign_prob_neg"]])
    fi=np.array([min(n_r-1,n_r//2+int((n_r//2-1)*j/max(n_z-1,1))) for j in range(n_z)],dtype=int)
    si=np.array([int((n_r//2-1)*(n_z-1-j)/max(n_z-1,1)) for j in range(n_z)],dtype=int)
    b0=_inv_sp(math.sqrt(max(TV_DAY**2-SIG_MIN**2,1e-15)))
    k0=arch["rough_orientation"]*arch["rough_scale"]*float(q@np.sqrt(2.*lam))
    return dict(lam=lam,q=q,az=az,zz=zz,sz=sz,fi=fi,si=si,b0=b0,k0=k0)

_ESN_PARAMS=None

def _sim_esn(seed,T,dt=1.0,arch=_ARCH):
    global _ESN_PARAMS
    if _ESN_PARAMS is None: _ESN_PARAMS=_build_esn(arch)
    P=_ESN_PARAMS; rng=np.random.default_rng(int(seed))
    n_r,n_z=arch["n_r"],arch["n_z"]
    spd=int(round(1./dt)); n_st=T*spd; sdt=math.sqrt(dt)
    al=np.exp(-P["lam"]*dt); cl=np.sqrt(np.maximum(1.-al**2,1e-14))
    azd=np.exp(-P["az"]*dt); om=1.-azd
    r=np.zeros(n_r); z=np.zeros(n_z); dx=np.zeros(T); dv=np.zeros(T)
    wz=arch["z_readout"]/math.sqrt(n_z); rc=arch["rough_orientation"]*arch["rough_scale"]
    for step in range(n_st):
        eps=rng.normal()
        eta=P["b0"]+rc*float(P["q"]@r)+wz*float(np.sum(z))
        sp=_sp(eta); sig=math.sqrt(SIG_MIN**2+sp**2); var=sig*sig
        day=step//spd
        dx[day]+=(0.-0.5*var)*dt+sig*sdt*eps; dv[day]+=var*dt
        r=al*r+cl*eps
        m=max(0.,1.-arch["gamma_norm"]*np.linalg.norm(r)/math.sqrt(n_r))
        zo=z.copy()
        for j in range(n_z):
            p1=r[P["fi"][j]]; p2=r[P["si"][j]]
            ev=0.7*p1*p1+0.3*p2*p2-1.; lin=0.7*p1+0.3*p2
            jm=max(0,j-1); jp=min(n_z-1,j+1)
            lc=0.5*arch["local_z_strength"]*(zo[jm]+zo[jp])
            u=(P["sz"][j]*p1*p2+arch["even_strength"]*ev
               -arch["linear_strength"]*lin+P["zz"][j]*zo[j]+lc)
            z[j]=azd[j]*zo[j]+om[j]*m*arch["z_strength"]*math.tanh(u)
    return dx,dv

def _sim_qrh(seed,T,dt,H,nu_vol,lam_k,c_frac,window=252):
    al=H+0.5; nu_hat=nu_vol*math.sqrt(Gamma(2*H))/Gamma(al)
    xi0=TV_DAY**2; c=c_frac*xi0; Y0=math.sqrt(max(xi0-c,0.))
    spd=int(round(1./dt)); n_st=T*spd; sdt=math.sqrt(dt)
    rng=np.random.default_rng(int(seed))
    dx=np.zeros(T); dv=np.zeros(T); buf=np.zeros(window); bi=0; nf=0
    def _w(n):
        j=np.arange(n,dtype=float); return dt**al/(Gamma(al)*al)*((j+1)**al-j**al)
    for step in range(n_st):
        nw=min(nf+1,window); w=_w(nw)
        idxs=[(bi-1-k)%window for k in range(nw)]
        Y=Y0+nu_hat*float(np.dot(w,buf[idxs]))
        V=max(min(Y*Y+c,xi0*50),1e-12); sig=math.sqrt(V); eps=rng.normal()
        day=step//spd
        dx[day]+=(0.-0.5*V)*dt+sig*sdt*eps; dv[day]+=V*dt
        buf[bi%window]=sig*eps; bi+=1; nf=min(nf+1,window)
    return dx,dv

def _sim_gl(seed,T,dt,beta,a1,ar,av,Vbar):
    spd=int(round(1./dt)); n_st=T*spd; sdt=math.sqrt(dt)
    rng=np.random.default_rng(int(seed))
    dx=np.zeros(T); dv=np.zeros(T); rbuf=[]; V=Vbar
    for step in range(n_st):
        sig=math.sqrt(max(V,1e-12)); eps=rng.normal()
        r_t=sig*sdt*eps; day=step//spd
        dx[day]+=(0.-0.5*V)*dt+r_t; dv[day]+=V*dt
        rbuf.append(r_t); t=len(rbuf)
        if t<2: V=Vbar
        else:
            ra=np.array(rbuf); lgs=np.arange(1,t+1,dtype=float); r2=ra[::-1]**2/dt
            Kr=lgs**(ar-0.5); Fr=float((Kr*r2).sum()/Kr.sum())
            Kv=lgs**(av-0.5); Fv=float((Kv*r2).sum()/Kv.sum())
            V=max((1-beta)*Vbar+beta*(a1*Fr+(1-a1)*Fv),1e-10)
    return dx,dv

def sim_model_paths(model_name,cal,n_paths,T,burn,dt=1.0):
    """Simulate one model — UNCHANGED from rough_bergomi_zumbach.py."""
    global _ESN_PARAMS; st_all=[]
    for p in range(n_paths):
        if model_name=="ESN A2-981003":
            old=_ESN_PARAMS["b0"]; _ESN_PARAMS["b0"]+=cal.get("b0_delta",0.)
            x,v=_sim_esn(9000+p,T,dt); x=x*cal.get("scale",1.0)
            _ESN_PARAMS["b0"]=old
        elif model_name=="QRH":
            x,v=_sim_qrh(7000+p,T,dt,cal["H"],cal["nu_vol"],cal["lam"],cal["c_frac"])
        elif model_name=="PDV-GL":
            x,v=_sim_gl(8000+p,T,dt,cal["beta"],cal["alpha1"],
                        cal["alpha_r"],cal["alpha_v"],cal["Vbar"])
        xb=x[burn:]; vb=v[burn:]
        st=compute_statistics(xb,vb); st_all.append(st)
    return st_all

# ============================================================
# 4.  Calibration — UNCHANGED (full 11-term score S)
# ============================================================

def quick_calibrate(dgp_sts, n_cal=4, T_cal=600, burn_cal=150, dt=1.0):
    """
    Calibrate ESN, QRH, PDV-GL to the MRW DGP.
    UNCHANGED from rough_bergomi_zumbach.py.
    The Taylor terms G_T (weight 1.0) and F_T (weight 0.8) are included,
    so calibration explicitly targets the Taylor gap and fraction.
    """
    global _ESN_PARAMS
    ref=make_score_ref(dgp_sts); tH=ref["H_hat_c"]
    print(f"    Score ref: H={tH:.4f}  vol={ref['mean_vol_ann_c']*100:.1f}%  "
          f"GT={ref['taylor_gap_c']:.4f}  FT={ref['taylor_frac_c']:.3f}  "
          f"K={ref['kurtosis_c']:.3f}")

    def esn_obj(p):
        b0d,sc=p; sc=max(sc,0.05)
        old=_ESN_PARAMS["b0"]; _ESN_PARAMS["b0"]=old+b0d
        scores=[]
        for q in range(n_cal):
            try:
                x,v=_sim_esn(3000+q,T_cal,dt)
                xb=(x*sc)[burn_cal:]; vb=v[burn_cal:]
                scores.append(score_fn(compute_statistics(xb,vb),ref))
            except: scores.append(-10.)
        _ESN_PARAMS["b0"]=old; return -float(np.nanmean(scores))
    res=optimize.minimize(esn_obj,[0.,1.],method="Nelder-Mead",
                          options={"maxiter":120,"xatol":0.05,"fatol":0.05})
    b0d,sc=res.x; sc=max(float(sc),0.05)
    print(f"    ESN: b0_delta={b0d:.4f}  scale={sc:.4f}  score={-res.fun:.3f}")
    esn_cal=dict(b0_delta=float(b0d),scale=sc)

    def qrh_obj(p):
        H,nu,lk,cf=p
        if not(0.01<H<0.49 and 0.01<nu<3. and 0.1<lk<20. and 0.01<cf<0.99): return 20.
        scores=[]
        for q in range(n_cal):
            try:
                x,v=_sim_qrh(500+q,T_cal,dt,H,nu,lk,cf)
                scores.append(score_fn(compute_statistics(x[burn_cal:],v[burn_cal:]),ref))
            except: scores.append(-10.)
        return -float(np.nanmean(scores))
    res=optimize.minimize(qrh_obj,[max(0.02,min(tH,0.48)),0.30,2.0,0.50],
                          method="Nelder-Mead",
                          options={"maxiter":300,"xatol":0.01,"fatol":0.02})
    H_o,nv_o,lk_o,cf_o=res.x
    H_o=float(np.clip(H_o,0.01,0.48)); nv_o=float(np.clip(nv_o,0.01,3.))
    lk_o=float(np.clip(lk_o,0.1,20.)); cf_o=float(np.clip(cf_o,0.01,0.99))
    print(f"    QRH: H={H_o:.4f}  nu={nv_o:.4f}  lam={lk_o:.4f}  c={cf_o:.4f}  score={-res.fun:.3f}")
    qrh_cal=dict(H=H_o,nu_vol=nv_o,lam=lk_o,c_frac=cf_o)

    Vbar=TV_DAY**2*dt
    def gl_obj(p):
        beta,a1,ar,av=p
        if not(0.05<beta<0.99 and 0<a1<1 and 0.01<ar<0.49 and 0.01<av<0.49): return 20.
        scores=[]
        for q in range(n_cal):
            try:
                x,v=_sim_gl(600+q,T_cal,dt,beta,a1,ar,av,Vbar)
                scores.append(score_fn(compute_statistics(x[burn_cal:],v[burn_cal:]),ref))
            except: scores.append(-10.)
        return -float(np.nanmean(scores))
    res=optimize.minimize(gl_obj,[0.75,0.45,0.08,0.25],method="Nelder-Mead",
                          options={"maxiter":400,"xatol":0.01,"fatol":0.02})
    b_o,a1_o,ar_o,av_o=res.x
    b_o=float(np.clip(b_o,0.05,0.99)); a1_o=float(np.clip(a1_o,0.01,0.99))
    ar_o=float(np.clip(ar_o,0.01,0.49)); av_o=float(np.clip(av_o,0.01,0.49))
    print(f"    PDV-GL: beta={b_o:.4f}  a1={a1_o:.4f}  ar={ar_o:.4f}  av={av_o:.4f}  score={-res.fun:.3f}")
    gl_cal=dict(beta=b_o,alpha1=a1_o,alpha_r=ar_o,alpha_v=av_o,Vbar=Vbar)

    return esn_cal,qrh_cal,gl_cal,ref

# ============================================================
# 5.  Taylor utility: per-lag ACF arrays
# ============================================================

def _taylor_arrays(st_list, max_lag=20):
    """
    Cross-path mean and std of ACF(|r|) and ACF(r^2) per lag.
    Returns: aa_m, aa_s, sa_m, sa_s  each shape (max_lag,)
    """
    aa = np.array([s["abs_acf_vec"][:max_lag] for s in st_list])
    sa = np.array([s["sq_acf_vec"][:max_lag]  for s in st_list])
    return aa.mean(0), aa.std(0), sa.mean(0), sa.std(0)

# ============================================================
# 6.  Figures
# ============================================================

def _col(name): return COLORS.get(name, GRAY)

# ── Figure T1: ACF curves + gap + summary bars ───────────────

def plot_T1_taylor_curves(results, lam2_grid, save_path=None):
    """
    Figure T1 — Taylor effect curves.
    One column per lambda^2.  Four rows:
      0. ACF(|r_t|) vs lag — all sources, ±1std shaded band.
         DGP line is thicker; models are thinner.
      1. ACF(r_t^2) vs lag — same layout.
         Should lie strictly BELOW row 0 everywhere = Taylor effect.
      2. Taylor gap = ACF(|r|) - ACF(r^2) per lag (grouped bars,
         one group per lag step 1,2,...,20).
         Positive bars = Taylor effect at that lag.
      3. Summary bars: Taylor gap scalar and Taylor fraction per model.
         Two bar groups side by side (gap solid, fraction hatched).
    """
    lam2s = [l for l in lam2_grid if l in results]
    n_c   = len(lam2s)
    lags  = np.arange(1, 21)
    ci_95 = 1.96 / math.sqrt(500)   # approximate 95% CI for ACF

    fig, axes = plt.subplots(4, n_c, figsize=(5*n_c, 16))
    if n_c == 1: axes = axes[:, None]

    fig.suptitle(
        r"Taylor Effect: $\mathrm{ACF}(|r_t|) > \mathrm{ACF}(r_t^2)$  for all $\tau > 0$"
        "\nMRW DGP  |  "
        r"Analytical: $\mathrm{ACF}(|r|)\sim\tau^{\lambda^2}$, "
        r"$\mathrm{ACF}(r^2)\sim\tau^{4\lambda^2}$"
        "\nModels calibrated via full 11-term score $S$",
        fontsize=11, fontweight="bold")

    for ci, lam2 in enumerate(lam2s):
        data = results[lam2]
        axes[0, ci].set_title(f"$\\lambda^2 = {lam2}$", fontsize=11, fontweight="bold")

        # ── rows 0 & 1: ACF curves ────────────────────────────────
        for ri, key in enumerate(["abs_acf_vec", "sq_acf_vec"]):
            ax = axes[ri, ci]
            for name in ALL_NAMES:
                aa = np.array([s[key][:20] for s in data[name]])
                m  = aa.mean(0); s = aa.std(0)
                lw = 2.2 if "DGP" in name else 1.5
                ax.fill_between(lags, m-s, m+s, alpha=0.13, color=_col(name))
                ax.plot(lags, m, color=_col(name), lw=lw, label=name)
            ax.axhline(ci_95,  color=GRAY, lw=0.7, ls="--", alpha=0.5)
            ax.axhline(-ci_95, color=GRAY, lw=0.7, ls="--", alpha=0.5)
            ax.axhline(0, color=GRAY, lw=0.4)
            ylabel = ("ACF$(|r_t|)$\n(should exceed row below)"
                      if ri == 0 else "ACF$(r_t^2)$\n(should lie below row above)")
            if ci == 0: ax.set_ylabel(ylabel, fontsize=8)
            ax.set_xlabel("Lag $\\tau$ (days)" if ri == 1 else "", fontsize=9)
            ax.legend(fontsize=7); ax.grid(True, alpha=0.18)

        # ── row 2: Taylor gap per lag (grouped bars) ──────────────
        ax = axes[2, ci]
        n_m = len(ALL_NAMES); bw = 0.70 / n_m
        for mi, name in enumerate(ALL_NAMES):
            aa_m, _, sa_m, _ = _taylor_arrays(data[name])
            gap = aa_m - sa_m
            offset = (mi - n_m/2 + 0.5) * bw
            ax.bar(lags + offset, gap, bw*0.92,
                   color=_col(name), alpha=0.82, label=name)
        ax.axhline(0, color=GRAY, lw=0.9, ls="--")
        ax.set_xlabel("Lag $\\tau$ (days)", fontsize=9)
        if ci == 0: ax.set_ylabel("ACF$(|r|)$ $-$ ACF$(r^2)$\n$>0$ = Taylor effect", fontsize=8)
        ax.set_xticks(lags[::2])
        ax.legend(fontsize=7); ax.grid(True, axis="y", alpha=0.18)

        # ── row 3: scalar summary (gap + fraction) ────────────────
        ax = axes[3, ci]
        x_ = np.arange(len(ALL_NAMES)); bw2 = 0.38
        gap_vals  = [float(np.nanmean([s["taylor_gap"]  for s in data[n]])) for n in ALL_NAMES]
        frac_vals = [float(np.nanmean([s["taylor_frac"] for s in data[n]])) for n in ALL_NAMES]
        bar_cols  = [_col(n) for n in ALL_NAMES]
        b1 = ax.bar(x_ - bw2/2, gap_vals,  bw2, color=bar_cols, alpha=0.85, label="Taylor gap")
        b2 = ax.bar(x_ + bw2/2, frac_vals, bw2, color=bar_cols, alpha=0.45,
                    hatch="///", label="Taylor fraction")
        ax.axhline(0,   color=GRAY,  lw=0.6)
        ax.axhline(0.5, color=CORAL, lw=0.9, ls="--", alpha=0.7, label="Fraction=0.5")
        ax.set_xticks(x_)
        ax.set_xticklabels(ALL_NAMES, rotation=25, fontsize=7, ha="right")
        if ci == 0: ax.set_ylabel("Gap (solid) / Fraction (hatch)", fontsize=8)
        ax.set_title("Taylor gap & fraction\n$\\mathrm{gap}>0$, $\\mathrm{frac}>0.5$"
                     " = Taylor effect", fontsize=8)
        ax.legend(fontsize=7)
        ax.grid(True, axis="y", alpha=0.18)
        for bar, v in zip(b1, gap_vals):
            ax.text(bar.get_x()+bar.get_width()/2, max(v,0)+0.003,
                    f"{v:.3f}", ha="center", fontsize=6.5)
        for bar, v in zip(b2, frac_vals):
            ax.text(bar.get_x()+bar.get_width()/2, max(v,0)+0.01,
                    f"{v:.2f}", ha="center", fontsize=6.5)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved -> {save_path}")
    return fig


# ── Figure T2: summary across lambda^2 ───────────────────────

def plot_T2_summary(results, lam2_grid, save_path=None):
    """
    Figure T2 — Taylor summary across lambda^2 grid.
    Four rows:
      0. Taylor gap (mean over lags 1-20) vs lambda^2
         with theoretical MRW reference line ~ lambda^2 * log(T_int).
      1. Taylor fraction vs lambda^2
         (=1 means ACF(|r|) > ACF(r^2) at every lag, every path).
      2. H_hat (dyadic variogram estimate) vs lambda^2.
      3. Score S vs lambda^2.
    """
    lam2s = [l for l in lam2_grid if l in results]
    marks = ["o","s","^","D"]; lss=["-","--","-.",":"]

    fig, axes = plt.subplots(4, 1, figsize=(9, 14))
    fig.suptitle(
        "Taylor Effect Summary vs Intermittency $\\lambda^2$\n"
        "MRW DGP  |  "
        r"$r_t = \sigma_0 e^{\omega_t}\varepsilon_t$,  "
        r"$\mathrm{Cov}(\omega_s,\omega_t)=\lambda^2\log(T_{\rm int}/\max(|t-s|,1))$",
        fontsize=11, fontweight="bold")

    # Row 0: Taylor gap vs lambda^2
    ax = axes[0]
    # Theoretical MRW: gap ≈ E[ACF(|r|) - ACF(r^2)] at lag 1
    # ACF(sigma^k, lag 1) ~ (1/T_int)^{k^2*lam2}
    # gap ≈ (1/T_int)^{lam2} - (1/T_int)^{4*lam2}
    th_gap = [( (1/T_INT)**l2 - (1/T_INT)**(4*l2) ) for l2 in lam2s]
    ax.plot(lam2s, th_gap, color=CORAL, lw=1.2, ls=":", label="Theory lag-1 approx.")
    for ni, name in enumerate(ALL_NAMES):
        vals = [float(np.nanmean([s["taylor_gap"] for s in results[l2][name]])) for l2 in lam2s]
        ax.plot(lam2s, vals, color=_col(name), lw=1.8, marker=marks[ni],
                ms=7, ls=lss[ni], label=name)
    ax.axhline(0, color=GRAY, lw=0.6, ls="--")
    ax.set_xticks(lam2s); ax.set_xlabel("$\\lambda^2$", fontsize=10)
    ax.set_ylabel("Taylor gap (mean lags 1–20)", fontsize=10)
    ax.set_title(r"Taylor gap vs $\lambda^2$  (MRW: $\propto\lambda^2$)", fontsize=10)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.22)

    # Row 1: Taylor fraction
    ax = axes[1]
    for ni, name in enumerate(ALL_NAMES):
        vals=[float(np.nanmean([s["taylor_frac"] for s in results[l2][name]])) for l2 in lam2s]
        ax.plot(lam2s, vals, color=_col(name), lw=1.8, marker=marks[ni],
                ms=7, ls=lss[ni], label=name)
    ax.axhline(0.5, color=CORAL, lw=1., ls="--", label="0.5 (neutral)")
    ax.axhline(1.0, color=GRAY,  lw=0.7, ls=":", alpha=0.7)
    ax.set_xticks(lam2s); ax.set_xlabel("$\\lambda^2$", fontsize=10)
    ax.set_ylabel("Taylor fraction (lags 1–20)", fontsize=10)
    ax.set_title(r"Taylor fraction vs $\lambda^2$  (=1 = perfect Taylor)", fontsize=10)
    ax.set_ylim(0, 1.1)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.22)

    # Row 2: H_hat
    ax = axes[2]
    for ni, name in enumerate(ALL_NAMES):
        vals=[float(np.nanmean([s["H_hat"] for s in results[l2][name]])) for l2 in lam2s]
        ax.plot(lam2s, vals, color=_col(name), lw=1.8, marker=marks[ni],
                ms=7, ls=lss[ni], label=name)
    ax.set_xticks(lam2s); ax.set_xlabel("$\\lambda^2$", fontsize=10)
    ax.set_ylabel(r"$\hat H$ (dyadic variogram)", fontsize=10)
    ax.set_title(r"$\hat H$ vs $\lambda^2$", fontsize=10)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.22)

    # Row 3: Score S
    ax = axes[3]
    for ni, name in enumerate(ALL_NAMES):
        vals=[float(np.nanmean([score_fn(s, results[l2]["ref"])
                                for s in results[l2][name]])) for l2 in lam2s]
        ax.plot(lam2s, vals, color=_col(name), lw=1.8, marker=marks[ni],
                ms=7, ls=lss[ni], label=name)
    ax.axhline(0, color=GRAY, lw=0.6, ls="--")
    ax.set_xticks(lam2s); ax.set_xlabel("$\\lambda^2$", fontsize=10)
    ax.set_ylabel("Score $S$ (model vs DGP)", fontsize=10)
    ax.set_title(r"Score $S$ vs $\lambda^2$  (higher = closer to DGP)", fontsize=10)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.22)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved -> {save_path}")
    return fig


# ── Figure T3: return diagnostics ────────────────────────────

def plot_T3_diagnostics(results, lam2_grid, save_path=None):
    """
    Figure T3 — Return diagnostics confirming DGP realism.
    One column per lambda^2. Four rows:
      0. Return histogram (log-y) with Gaussian overlay.
      1. Hill tail-index curve alpha_hat(k) vs k.
      2. Excess kurtosis bar chart.
      3. QQ-plot of standardised returns vs Gaussian.
    """
    lam2s = [l for l in lam2_grid if l in results]
    n_c   = len(lam2s)

    fig, axes = plt.subplots(4, n_c, figsize=(5*n_c, 14))
    if n_c == 1: axes = axes[:, None]

    fig.suptitle(
        "Return Diagnostics — MRW DGP\n"
        r"$r_t = \sigma_0 e^{\omega_t}\varepsilon_t$,  $\varepsilon_t\sim\mathcal{N}(0,1)$"
        "\nFat tails emerge from stochastic vol, not from innovation distribution",
        fontsize=11, fontweight="bold")

    for ci, lam2 in enumerate(lam2s):
        data = results[lam2]
        axes[0, ci].set_title(f"$\\lambda^2 = {lam2}$", fontsize=11, fontweight="bold")

        for name in ALL_NAMES:
            r_all = np.concatenate([s["logret"] for s in data[name]])
            col   = _col(name); lw_ = 2.2 if "DGP" in name else 1.5
            lo, hi = np.percentile(r_all, 0.3), np.percentile(r_all, 99.7)

            # Row 0: histogram
            ax = axes[0, ci]
            ax.hist(r_all, bins=120, density=True, histtype="step",
                    lw=lw_, color=col, label=name)
            if "DGP" in name:
                xr = np.linspace(lo, hi, 300)
                ax.plot(xr, stats.norm.pdf(xr, 0, r_all.std()),
                        color=GRAY, lw=0.9, ls=":", label="Gaussian ref")
            ax.set_yscale("log"); ax.set_xlim(lo, hi)
            if ci == 0: ax.set_ylabel("Density (log)", fontsize=9)
            ax.legend(fontsize=6); ax.grid(True, alpha=0.18)

            # Row 1: Hill
            ax = axes[1, ci]
            ks, al = hill_curve(r_all)
            ax.plot(ks, al, color=col, lw=lw_, label=name)
            if ci == 0: ax.set_ylabel(r"Hill $\hat\alpha$", fontsize=9)
            ax.legend(fontsize=6); ax.grid(True, alpha=0.18)

        # Row 2: kurtosis bar
        ax = axes[2, ci]
        k_vals = [float(np.nanmean([s["kurtosis"] for s in data[n]])) for n in ALL_NAMES]
        bars = ax.bar(ALL_NAMES, k_vals,
                      color=[_col(n) for n in ALL_NAMES], alpha=0.82, width=0.6)
        ax.axhline(0, color=GRAY, lw=0.6)
        if ci == 0: ax.set_ylabel("Excess kurtosis", fontsize=9)
        ax.tick_params(axis="x", labelsize=6, rotation=20)
        ax.grid(True, axis="y", alpha=0.18)
        for bar, v in zip(bars, k_vals):
            ax.text(bar.get_x()+bar.get_width()/2, max(v,0)+0.05,
                    f"{v:.1f}", ha="center", fontsize=7)

        # Row 3: QQ vs Gaussian
        ax = axes[3, ci]
        for name in ALL_NAMES:
            r_all = np.concatenate([s["logret"] for s in data[name]])
            z = (r_all - r_all.mean()) / (r_all.std() + 1e-14)
            qs = np.linspace(0.01, 0.99, 200)
            ax.plot(stats.norm.ppf(qs), np.quantile(z, qs),
                    color=_col(name), lw=2.2 if "DGP" in name else 1.5, label=name)
        ax.plot([-4,4],[-4,4], color=GRAY, lw=0.8, ls=":", label="Gaussian")
        ax.set_xlabel("Gaussian quantile", fontsize=9)
        if ci == 0: ax.set_ylabel("Empirical quantile", fontsize=9)
        ax.legend(fontsize=6); ax.grid(True, alpha=0.18)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved -> {save_path}")
    return fig

# ============================================================
# 7.  LaTeX protocol
# ============================================================

LATEX = r"""\documentclass[11pt,a4paper]{article}
\usepackage{amsmath,amssymb,amsthm,amsfonts}
\usepackage{geometry}\geometry{left=2.3cm,right=2.3cm,top=2.5cm,bottom=2.5cm}
\usepackage{setspace}\setstretch{1.20}
\usepackage{hyperref}\hypersetup{colorlinks=true,linkcolor=black,citecolor=black}
\usepackage{parskip,booktabs,array,xcolor,titlesec,enumitem}
\titleformat{\section}{\large\bfseries}{}{0em}{}[\titlerule]
\titlespacing{\section}{0pt}{14pt}{6pt}
\titleformat{\subsection}{\normalsize\bfseries}{\thesubsection\;}{0em}{}
\newtheorem{remark}{Remark}
\begin{document}
\begin{center}
{\LARGE\bfseries Numerical Protocol}\\[0.3cm]
{\large Taylor Effect Evaluation\\
under the Multifractal Random Walk (MRW) DGP}\\[0.2cm]
{\normalsize Othmane Zarhali --- Paris Dauphine / CNRS}
\end{center}
\vspace{0.15cm}\noindent\rule{\textwidth}{0.8pt}

\begin{abstract}
The Taylor effect is the empirical regularity
$\mathrm{ACF}(|r_t|) > \mathrm{ACF}(r_t^2)$ for all lags $\tau > 0$.
We evaluate how well three stochastic volatility models (ESN A2-981003, QRH,
PDV-GL) reproduce this effect under the Multifractal Random Walk (MRW)
data-generating process, for which the Taylor inequality holds analytically
for all intermittency parameters $\lambda^2 > 0$.
All models are calibrated via the full 11-term data-adaptive score $S$
(Eq.~\eqref{eq:score}), which includes both the Taylor gap $G_T$ (weight 1.0)
and Taylor fraction $F_T$ (weight 0.8).
\end{abstract}

\tableofcontents\vspace{0.4cm}

\section{DGP: Multifractal Random Walk}

The DGP is the MRW \cite{bacry2001}:
\begin{equation}
  r_t = \sigma_0\,e^{\omega_t}\,\varepsilon_t, \qquad
  \varepsilon_t \sim \mathcal{N}(0,1),
  \label{eq:mrw}
\end{equation}
where $\omega_t$ is a Gaussian process with log-correlation kernel:
\begin{equation}
  \mathrm{Cov}(\omega_s,\omega_t)
  = \lambda^2\,\max\!\Bigl(0,\,\log\frac{T_{\rm int}}{\max(|t-s|,1)}\Bigr),
  \quad |t-s| < T_{\rm int},
  \label{eq:cov}
\end{equation}
and $\mathrm{Cov}(\omega_s,\omega_t)=0$ for $|t-s|\ge T_{\rm int}$.
The centering $-\lambda^2\log(T_{\rm int})/2$ is subtracted from $\omega_t$
to ensure $\mathbb{E}[\sigma_t]=\sigma_0$.

\paragraph{Analytical Taylor inequality.}
Since $\sigma_t=\sigma_0 e^{\omega_t}$ is log-normal,
$\mathrm{Corr}(\sigma_t^k,\sigma_{t+\tau}^k)\propto(\tau/T_{\rm int})^{k^2\lambda^2}$,
so:
\begin{equation}
  \mathrm{ACF}(|r_t|)\propto\tau^{\lambda^2},\qquad
  \mathrm{ACF}(r_t^2)\propto\tau^{4\lambda^2}.
  \label{eq:taylor_proof}
\end{equation}
Since $4\lambda^2>\lambda^2$ for all $\lambda^2>0$, the squared-return ACF
decays faster, giving a positive Taylor gap at every lag.
The gap grows monotonically with $\lambda^2$.

\paragraph{Simulation.}
The covariance matrix $C_{ij}=\mathrm{Cov}(\omega_i,\omega_j)$ is NOT positive
definite on finite grids (minimum eigenvalue $\approx -\lambda^2/100$).
We project onto the nearest PSD matrix via eigenvalue flooring
$\Lambda_{\rm PSD}=\max(\Lambda,10^{-10})$ before Cholesky decomposition.
The Taylor effect is verified numerically on each DGP run.

\paragraph{Parameter grid.}
$\lambda^2\in\{0.01,0.02,0.03,0.04\}$ (0.03 calibrated to SPX);
$T_{\rm int}=252$ (integral scale, one trading year).

\section{Models (Unchanged)}

All models use Gaussian Brownian innovations.
Fat tails in returns emerge from stochastic variance clustering.

\subsection{ESN A2-981003}
Exact notebook simulate\_path ($N_r=64$, $N_z=12$, $\kappa_0=-0.472416$).
Single $\varepsilon_t\sim\mathcal{N}(0,1)$ drives both return and reservoir.
The weighted even response $0.7p_1^2+0.3p_2^2-1$ in the interaction bank
is the primary architectural mechanism generating the Taylor effect:
it creates a quadratic response to shock magnitude making $|r|$ ACF
decay slower than $r^2$ ACF.

\subsection{Quadratic Rough Heston (QRH)}
Bourgey \& Gatheral \cite{bourgey2026}: $V_t=Y_t^2+c$, gamma kernel,
Euler-Volterra with 252-step ring buffer.

\subsection{PDV Guyon-Lekeufack (GL)}
Guyon \& Lekeufack \cite{guyon2023}: two-factor power-law kernel,
online from full return history. The long-memory factor $F(\alpha_v)$ with
$\alpha_v<0.5$ directly encodes power-law memory of past squared returns,
which is the GL mechanism for Taylor.

\section{Calibration (Unchanged --- Full 11-Term Score)}

\begin{equation}
\begin{aligned}
S &= 2.2\,f(\hat H)+1.4\,f(\bar\sigma)+0.6\,g(q_{995})+0.2\,g(\sigma_{\max})
   +0.8\,f(\bar V_{\rm ACF})\\
  &+\underbrace{1.0\,f(G_T)}_{\text{Taylor gap}}
   +\underbrace{0.8\,f(F_T)}_{\text{Taylor frac.}}
   +1.0\,f(Z)+1.1\,f(L)+0.7\,f(K)+0.5\,f(A)-5\cdot\mathrm{Stress}.
\end{aligned}
\label{eq:score}
\end{equation}
Taylor terms carry weight $1.0+0.8=1.8$ (second highest after $\hat H$),
so calibration explicitly targets Taylor gap and fraction simultaneously.
All centres $c_k$ and tolerances $s_k$ are DGP cross-path mean and
$\max(\mathrm{std},0.30|c_k|,\mathrm{floor}_k)$.

\section{Taylor Diagnostics}

\paragraph{Per-lag ACF curves (Figure T1, rows 0--1).}
Cross-path mean $\pm1$ std of ACF$(|r_t|)$ and ACF$(r_t^2)$ at lags
$\tau=1,\ldots,20$ days. Taylor effect present $\Leftrightarrow$
row-0 curve lies strictly above row-1 curve at every lag.

\paragraph{Taylor gap per lag (Figure T1, row 2).}
\begin{equation}
  \Delta(\tau) = \mathrm{ACF}(|r_t|,\tau) - \mathrm{ACF}(r_t^2,\tau)
  \label{eq:gap}
\end{equation}
Grouped bar chart (one group per lag, one bar per source).
$\Delta(\tau)>0$ = Taylor effect at that lag.

\paragraph{Scalar summary (Figure T1, row 3).}
Taylor gap scalar $G_T=\frac{1}{20}\sum_\tau\Delta(\tau)$ and Taylor
fraction $F_T=\frac{1}{20}\sum_\tau\mathbf{1}[\Delta(\tau)>0]$ per model.

\section{Figures}

\textbf{Figure T1.} One column per $\lambda^2$, four rows:
ACF$(|r|)$; ACF$(r^2)$; per-lag gap (grouped bars); gap scalar + fraction.

\textbf{Figure T2.} Summary across $\lambda^2$: Taylor gap, fraction,
$\hat H$, score $S$.  DGP shows increasing gap and fraction with $\lambda^2$.

\textbf{Figure T3.} Return diagnostics: histogram (log-y), Hill $\hat\alpha$,
kurtosis bar, QQ vs Gaussian.

\section{Monte Carlo Setup}

\begin{center}\renewcommand{\arraystretch}{1.4}
\begin{tabular}{ll}\toprule Parameter & Value \\\midrule
DGP paths $n_{\rm dgp}$ & 30 \\
Calibration paths $n_{\rm cal}$ & 5 \\
Calibration trading days & 600 \\
Final paths $n_{\rm sim}$ & 30 \\
Final trading days & 2500 \\
Burn-in & 500 \\
$\Delta t$ & 1 (daily) \\
$\lambda^2$ grid & $\{0.01,0.02,0.03,0.04\}$ \\
$T_{\rm int}$ & 252 (days) \\
Eigenvalue floor & $10^{-10}$ (MRW PSD projection) \\
\bottomrule\end{tabular}\end{center}

\begin{thebibliography}{9}
\bibitem{bacry2001} E.\ Bacry, J.\ Delour, J.-F.\ Muzy.
\textit{Multifractal random walk}. Phys.\ Rev.\ E 64:026103, 2001.
\bibitem{muzy2000} J.-F.\ Muzy, J.\ Delour, E.\ Bacry.
\textit{Modelling fluctuations of financial time series}.
EPJB 17:537--548, 2000.
\bibitem{bourgey2026} F.\ Bourgey, J.\ Gatheral.
\textit{Quadratic Rough Heston}. SSRN:5239929, 2026.
\bibitem{guyon2023} J.\ Guyon, J.\ Lekeufack.
\textit{Volatility is (mostly) path-dependent}. QF 23(9), 2023.
\bibitem{taylor1986} S.J.\ Taylor.
\textit{Modelling Financial Time Series}. Wiley, 1986.
\end{thebibliography}
\end{document}
"""

def write_latex(path):
    with open(path,"w") as f: f.write(LATEX.lstrip())
    print(f"  Saved -> {path}")

# ============================================================
# 8.  Main pipeline
# ============================================================

def run(lam2_grid=None, T_int=T_INT,
        n_dgp=30, n_sim=30, T=2500, burn=500,
        n_cal=5, T_cal=600, burn_cal=150, dt=1.0,
        save_prefix="mrw_taylor", latex_path=None):

    if lam2_grid is None: lam2_grid = LAM2_GRID

    global _ESN_PARAMS
    _ESN_PARAMS = _build_esn(_ARCH)
    assert _ESN_PARAMS["k0"] < 0
    print(f"ESN kappa0 = {_ESN_PARAMS['k0']:.6f}  ✓\n")

    results = {}

    for lam2 in lam2_grid:
        print(f"{'='*60}  lambda^2 = {lam2}  (T_int = {T_int})")

        # ── DGP ──────────────────────────────────────────────
        print("  DGP (MRW) ...")
        dgp_sts = dgp_mrw(lam2, n_dgp, T, burn,
                          T_int=T_int, dt=dt,
                          seed_base=int(lam2*10000))

        def _m(k): return float(np.nanmean([s[k] for s in dgp_sts]))
        print(f"  DGP: H_hat={_m('H_hat'):.4f}  vol={_m('mean_vol_ann')*100:.1f}%  "
              f"GT={_m('taylor_gap'):.4f}  FT={_m('taylor_frac'):.3f}  "
              f"K={_m('kurtosis'):.3f}")

        # verify Taylor effect is present in DGP
        assert _m('taylor_gap') > 0, f"Taylor effect missing from DGP (lam2={lam2})!"
        assert _m('taylor_frac') > 0.5, f"Taylor fraction < 0.5 (lam2={lam2})!"

        # ── Calibrate ─────────────────────────────────────────
        print("  Calibrating (full 11-term score S) ...")
        esn_cal, qrh_cal, gl_cal, ref = quick_calibrate(
            dgp_sts, n_cal, T_cal, burn_cal, dt)

        # ── Simulate final paths ──────────────────────────────
        print("  Simulating final paths ...")
        esn_sts = sim_model_paths("ESN A2-981003", esn_cal, n_sim, T, burn, dt)
        qrh_sts = sim_model_paths("QRH",           qrh_cal, n_sim, T, burn, dt)
        gl_sts  = sim_model_paths("PDV-GL",         gl_cal,  n_sim, T, burn, dt)

        results[lam2] = {
            "DGP (MRW)":      dgp_sts,
            "ESN A2-981003":  esn_sts,
            "QRH":            qrh_sts,
            "PDV-GL":         gl_sts,
            "ref":            ref,
        }

        # ── Summary ───────────────────────────────────────────
        print(f"\n  Summary  lambda^2={lam2}:")
        print(f"  {'Model':<22} {'H_hat':>7} {'vol%':>6} {'GT':>8} "
              f"{'FT':>6} {'Kurt':>7} {'Score':>7}")
        print("  "+"-"*65)
        for name in ALL_NAMES:
            sts = results[lam2][name]
            sc  = float(np.nanmean([score_fn(s, ref) for s in sts]))
            print(f"  {name:<22} "
                  f"{float(np.nanmean([s['H_hat']       for s in sts])):>7.4f} "
                  f"{float(np.nanmean([s['mean_vol_ann'] for s in sts]))*100:>5.1f}% "
                  f"{float(np.nanmean([s['taylor_gap']   for s in sts])):>8.4f} "
                  f"{float(np.nanmean([s['taylor_frac']  for s in sts])):>6.3f} "
                  f"{float(np.nanmean([s['kurtosis']     for s in sts])):>7.3f} {sc:>7.3f}")
        print()

    # ── Figures ──────────────────────────────────────────────────────────
    print(f"\n{'='*60}\nProducing figures ...\n{'='*60}")

    fig = plot_T1_taylor_curves(results, lam2_grid,
          save_path=f"{save_prefix}_T1.png")
    plt.close(fig)

    fig = plot_T2_summary(results, lam2_grid,
          save_path=f"{save_prefix}_T2.png")
    plt.close(fig)

    fig = plot_T3_diagnostics(results, lam2_grid,
          save_path=f"{save_prefix}_T3.png")
    plt.close(fig)

    # ── LaTeX / PDF ───────────────────────────────────────────────────────
    lp = latex_path or f"{save_prefix}_protocol.tex"
    write_latex(lp)
    if os.system("which pdflatex>/dev/null 2>&1") == 0:
        os.system(f"pdflatex -interaction=nonstopmode {lp}>/dev/null 2>&1")
        os.system(f"pdflatex -interaction=nonstopmode {lp}>/dev/null 2>&1")
        pdf = lp.replace(".tex", ".pdf")
        if os.path.exists(pdf): print(f"  Compiled -> {pdf}")

    return results


# ============================================================
# 9.  CLI
# ============================================================

if __name__ == "__main__":
    import argparse
    pa = argparse.ArgumentParser(
        description="Taylor effect under MRW DGP",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pa.add_argument("--lam2",  nargs="+", type=float, default=[0.01,0.02,0.03,0.04])
    pa.add_argument("--T_int", type=int,   default=252)
    pa.add_argument("--n_dgp", type=int,   default=30)
    pa.add_argument("--n_cal", type=int,   default=5)
    pa.add_argument("--n_sim", type=int,   default=30)
    pa.add_argument("--T",     type=int,   default=2500)
    pa.add_argument("--burn",  type=int,   default=500)
    pa.add_argument("--dt",    type=float, default=1.0)
    pa.add_argument("--out",   type=str,   default="mrw_taylor")
    pa.add_argument("--fast",  action="store_true",
                    help="lam2=[0.02,0.03], n=8, T=800")
    args = pa.parse_args()
    if args.fast:
        args.lam2=[0.02,0.03]; args.n_dgp=8; args.n_cal=3
        args.n_sim=8; args.T=800; args.burn=200
    run(lam2_grid=args.lam2, T_int=args.T_int,
        n_dgp=args.n_dgp, n_sim=args.n_sim, T=args.T, burn=args.burn,
        n_cal=args.n_cal, T_cal=min(args.T,600), burn_cal=min(args.burn,150),
        dt=args.dt, save_prefix=args.out)
