"""
roughness_evaluation.py
=======================
DGP: log-normal SV with fractional Brownian motion.
    r_t = sigma_t * eps_t,   eps_t ~ N(0,1)
    sigma_t = sigma_0 * exp( lambda * B^H(t) - lambda^2/2 * t^{2H} )
Models (UNCHANGED): ESN A2-981003, QRH, PDV-GL — all Gaussian.
Calibration (UNCHANGED): score S_H = f(H_hat_model, H_hat_DGP, s_H).
Diagnostics: H_hat estimation + CI, ACov vs tau^{2H_hat}, scale invariance.
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

TRADING_DAYS = 252.0
TV_ANN  = 0.20
TV_DAY  = TV_ANN / math.sqrt(TRADING_DAYS)
SIG_MIN = 0.01 / math.sqrt(TRADING_DAYS)
H_GRID      = [0.05, 0.10, 0.40]   # reduced to 3 values per request
LAMBDA_GRID = [0.50, 1.00, 2.00]
Q_GRID      = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
DYADIC_LAGS = [1, 2, 4, 8, 16, 32, 64]
COLORS = {"DGP (fBM SV)":"#888780","ESN A2-981003":"#1a4d80","QRH":"#D85A30","PDV-GL":"#1D9E75"}
CORAL="#c0392b"; GRAY="#888780"
MODEL_NAMES = ["ESN A2-981003","QRH","PDV-GL"]
_ARCH = dict(matrix_seed=202695547565,n_r=96,n_z=12,H_target=0.08,
             z_strength=0.19288535085059877,
             even_strength=2.2291053209265703,
             linear_strength=0.1930215632408165, gamma_norm=0.9651266728165394,
             local_z_strength=0.04759055003377607, zz_scale=0.039730814766135034,
             sign_prob_neg=0.22189311131648579, rough_orientation=-1.0)
# _ARCH is the ESN model's current optimal architecture: Nr=96, Nz=12,
# with 7 shared hyperparameters. H_target is used only as a start-up
# default; H, lam_lo, lam_hi, az_lo, az_hi, rough_scale, zr_lo, zr_hi,
# m1, m2 are all CALIBRATED per (H,lambda) DGP cell below
# (quick_calibrate), not fixed. The read-out is
#   eta_t = b0 + rho_r*rough_scale*(q^T r_t) + sum_j w_{j,z} z_{j,t}
#           + r_t^T Q r_t,
#   w_{j,z} = z_readout_j/sqrt(Nz), z_readout_j = geomspace(zr_lo,zr_hi,Nz)[j],
#   Q = (1/Nr)*(m1*I + m2*q q^T)  (symmetric, no PD constraint).
ROUGH_SCALE_BOUNDS = (0.05, 2.0)
ZR_LO_BOUNDS = (0.01, 0.5)
ZR_HI_BOUNDS = (0.01, 0.5)
M1_BOUNDS = (-2.0, 2.0)
M2_BOUNDS = (-2.0, 2.0)

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

# ── Full stylized-fact statistics (same as student_t_fat_tail_v3.py) ───────

def compute_statistics(daily_x, daily_var):
    """
    Eleven stylized-fact statistics used by the global score S.
    Mirrors compute_statistics() in student_t_fat_tail_v3.py exactly.
    """
    x = np.asarray(daily_x, float)
    v = np.asarray(daily_var, float)
    sigma_daily = np.sqrt(np.maximum(v, 1e-30))
    lv_         = np.log(v + 1e-30)
    lags        = np.arange(1, 21)

    ret_acf = np.array([_corr(x[:-L], x[L:]) for L in lags])
    abs_acf = np.array([_corr(np.abs(x[:-L]), np.abs(x[L:])) for L in lags])
    sq_acf  = np.array([_corr(x[:-L]**2, x[L:]**2) for L in lags])
    vol_acf = np.array([_corr(v[:-L], v[L:]) for L in lags])
    lev     = np.array([_corr(x[:-L], v[L:]) for L in lags])

    xs, ys = [], []
    for L in [1,2,4,8,16,32,64]:
        if L < len(lv_)//4:
            d = lv_[L:]-lv_[:-L]; vv = float(np.mean(d*d))
            if vv > 1e-30 and np.isfinite(vv):
                xs.append(math.log(L)); ys.append(math.log(vv))
    H = 0.5*float(np.polyfit(xs,ys,1)[0]) if len(xs)>=3 else np.nan

    zvals = []
    for L in [5,10,20]:
        n = len(x)-2*L
        if n > 30:
            csx = np.concatenate([[0.], np.cumsum(x)])
            csv = np.concatenate([[0.], np.cumsum(v)])
            idx = L+np.arange(n)
            pR = csx[idx]-csx[idx-L]; fR = csx[idx+L]-csx[idx]
            pV = (csv[idx]-csv[idx-L])/L; fV = (csv[idx+L]-csv[idx])/L
            zvals.append(_corr(pR**2,fV)-_corr(pV,fR**2))

    c = x-x.mean(); var = float(np.var(c))
    kurt = float(np.mean(c**4)/(var**2+1e-30)-3.)
    ann  = math.sqrt(TRADING_DAYS)
    return dict(
        H_hat        = float(H),
        mean_vol_ann = float(sigma_daily.mean())*ann,
        q995_vol_ann = float(np.quantile(sigma_daily,0.995))*ann,
        max_vol_ann  = float(sigma_daily.max())*ann,
        mean_vol_acf = float(np.nanmean(vol_acf)),
        taylor_gap   = float(np.nanmean(abs_acf)-np.nanmean(sq_acf)),
        taylor_frac  = float(np.nanmean(abs_acf>sq_acf)),
        zumbach      = float(np.nanmean(zvals)),
        leverage     = float(np.nanmean(lev)),
        kurtosis     = kurt,
        max_ret_acf  = float(np.nanmax(np.abs(ret_acf))),
    )

def make_score_ref(dgp_stats):
    """
    Data-adaptive score reference: c_k = DGP mean, s_k = max(std, 0.30|c|, floor).
    Identical to make_score_ref() in student_t_fat_tail_v3.py.
    """
    def _m(k):  return float(np.nanmean([s[k] for s in dgp_stats]))
    def _sd(k): return float(np.nanstd( [s[k] for s in dgp_stats]))
    floors = dict(H_hat=0.01, mean_vol_ann=0.01, q995_vol_ann=0.05,
                  max_vol_ann=0.10, mean_vol_acf=0.005, taylor_gap=0.002,
                  taylor_frac=0.05, zumbach=0.005, leverage=0.005,
                  kurtosis=0.20, max_ret_acf=0.005)
    ref = {}
    for k, fl in floors.items():
        c = _m(k); s = _sd(k)
        ref[k+"_c"] = c
        ref[k+"_s"] = max(s, abs(c)*0.30, fl)
    ref["stress_mx"] = max(ref["max_vol_ann_c"]*2.0, 1.5)
    return ref

SCORE_MAX = 2.2+1.4+0.6+0.2+0.8+1.0+0.8+1.0+1.1+0.7+0.5 + 5.0   # = 15.3, worst possible S

def score_fn(st, ref):
    """
    Full data-adaptive score S (11 terms), SMOOTH (Gaussian-kernel)
    version -- EXACT convention of the ESN hyperparameter-search
    protocol (V5): S >= 0, MINIMISED, S = 0 <=> a perfect match on every
    one of the 11 stylised facts and no stress event.
    S = sum_k w_k*(1-f_k) + 5*Stress, with
        f_smooth(x,c,s) = exp(-0.5*((x-c)/s)^2)             in (0,1]
        g_smooth(x,c,s) = exp(-0.5*(max(0,x-c)/s)^2)         in (0,1]
    replacing the old piecewise-linear tent max(0,1-|x-c|/s), which hits
    EXACTLY 0 (zero gradient) beyond one tolerance width -- the smooth
    kernel never flattens, giving Nelder-Mead a usable gradient at any
    distance from target (see ESN hyperparameter-search protocol V5, the
    "Smooth Score Function" section, for the full derivation).
    """
    def f(x,c,s): return math.exp(-0.5*((x-c)/max(s,1e-8))**2)
    def g(x,c,s): return math.exp(-0.5*(max(0.,x-c)/max(s,1e-8))**2)
    H   = st["H_hat"];          vol = st["mean_vol_ann"]
    q995= st["q995_vol_ann"];   mx  = st["max_vol_ann"]
    V   = st["mean_vol_acf"];   GT  = st["taylor_gap"]
    FT  = st["taylor_frac"];    Z   = st["zumbach"]
    L   = st["leverage"];       K   = st["kurtosis"]
    A   = st["max_ret_acf"]
    stress = int(mx > ref["stress_mx"] or vol < 0.05 or vol > 1.50)
    d  = 2.2*(1.-f(H,   ref["H_hat_c"],       ref["H_hat_s"]))
    d += 1.4*(1.-f(vol, ref["mean_vol_ann_c"], ref["mean_vol_ann_s"]))
    d += 0.6*(1.-g(q995,ref["q995_vol_ann_c"], ref["q995_vol_ann_s"]))
    d += 0.2*(1.-g(mx,  ref["max_vol_ann_c"],  ref["max_vol_ann_s"]))
    d += 0.8*(1.-f(V,   ref["mean_vol_acf_c"], ref["mean_vol_acf_s"]))
    d += 1.0*(1.-f(GT,  ref["taylor_gap_c"],   ref["taylor_gap_s"]))
    d += 0.8*(1.-f(FT,  ref["taylor_frac_c"],  ref["taylor_frac_s"]))
    d += 1.0*(1.-f(Z,   ref["zumbach_c"],      ref["zumbach_s"]))
    d += 1.1*(1.-f(L,   ref["leverage_c"],     ref["leverage_s"]))
    d += 0.7*(1.-f(K,   ref["kurtosis_c"],     ref["kurtosis_s"]))
    d += 0.5*(1.-f(A,   ref["max_ret_acf_c"],  ref["max_ret_acf_s"]))
    d += 5.0*stress
    return float(d)

def _stat_from_path(x, v):
    """Thin wrapper: compute full score stats from (x,v) arrays."""
    return compute_statistics(x, v)

# ── Hurst estimation and roughness diagnostics (unchanged) ──────────────────

def hurst_ols(lv, lags=None):
    if lags is None: lags=DYADIC_LAGS
    xs,ys=[],[]
    for L in lags:
        if L<len(lv)//4:
            d=lv[L:]-lv[:-L]; vv=float(np.mean(d*d))
            if vv>1e-30 and np.isfinite(vv):
                xs.append(math.log(L)); ys.append(math.log(vv))
    if len(xs)<3: return np.nan,np.nan,np.nan,np.nan,np.nan
    xs,ys=np.array(xs),np.array(ys)
    sl,ic,r,_,se=stats.linregress(xs,ys)
    return sl/2.0,sl,ic,r**2,se/2.0

def autocovariance_vec(lv,max_lag):
    xc=lv-lv.mean(); C0=float(np.mean(xc**2))
    Ck=np.array([float(np.mean(xc[:-k]*xc[k:])) for k in range(1,max_lag+1)])
    return C0,Ck

def structure_function(lv,q_grid,lags=None):
    if lags is None: lags=DYADIC_LAGS
    SF=np.zeros((len(lags),len(q_grid)))
    for li,L in enumerate(lags):
        if L<len(lv)//4:
            d=np.abs(lv[L:]-lv[:-L])
            for qi,q in enumerate(q_grid):
                SF[li,qi]=float(np.mean(d**q))
        else: SF[li,:]=np.nan
    return SF

def generalised_hurst(lv,q_grid,lags=None):
    if lags is None: lags=DYADIC_LAGS
    SF=structure_function(lv,q_grid,lags)
    ll=np.log(np.array(lags,float))
    zeta=np.full(len(q_grid),np.nan)
    for qi in range(len(q_grid)):
        sf=SF[:,qi]; ok=np.isfinite(sf)&(sf>0)
        if ok.sum()>=3:
            sl,*_=stats.linregress(ll[ok],np.log(sf[ok]))
            zeta[qi]=float(sl)
    Hq=zeta/np.array(q_grid,float)
    z1=zeta[np.argmin(np.abs(np.array(q_grid)-1.0))]
    mu=zeta-np.array(q_grid,float)*z1
    return zeta,Hq,mu

def compute_roughness_stats(lv,max_lag=40,q_grid=None):
    if q_grid is None: q_grid=Q_GRID
    H_hat,sl,ic,r2,H_se=hurst_ols(lv)
    C0,Ck=autocovariance_vec(lv,max_lag)
    zeta,Hq,mu=generalised_hurst(lv,q_grid)
    return dict(H_hat=H_hat,H_se=H_se,r2_variogram=r2,
                C0=C0,Ck=Ck,max_lag=max_lag,zeta=zeta,Hq=Hq,mu=mu,lv=lv)

def acf_scaling_regression(C0,Ck,H_hat,short_lag=20):
    if C0<1e-20 or not np.isfinite(H_hat):
        return dict(a=np.nan,b=np.nan,r2=np.nan,tau_2H=np.zeros(len(Ck)),rho=np.zeros(len(Ck)),n_fit=0)
    rho=Ck/C0; n=len(rho)
    tau=np.arange(1,n+1,dtype=float); tau_2H=tau**(2.0*H_hat)
    nf=min(short_lag,n)
    X=np.column_stack([np.ones(nf),tau_2H[:nf]])
    coef,*_=np.linalg.lstsq(X,rho[:nf],rcond=None); a,b=coef
    yfit=X@coef; ss_res=((rho[:nf]-yfit)**2).sum()
    ss_tot=((rho[:nf]-rho[:nf].mean())**2).sum()
    r2=1.0-ss_res/max(ss_tot,1e-12)
    return dict(a=float(a),b=float(b),r2=float(r2),tau_2H=tau_2H,rho=rho,n_fit=nf)

# ── DGP ─────────────────────────────────────────────────────
# *** FIX: annualised time argument for the fBM (see ESN hyperparameter-
# search protocol V1/V5, "fBM-persist data-generation collapse") ***
# The fBM covariance and the Ito centring term -lambda^2/2 * t^{2H} must
# use t measured in the SAME units as lambda (an O(1) annualised
# vol-of-vol parameter), i.e. YEARS, not raw day-count. Using raw day
# count t=1..T (T up to 2500) makes the centring term grow like
# T^{2H}: mild for H=0.05/0.10 (T^{0.1}/T^{0.2} stays small) but
# catastrophic for H=0.40 (T^{0.8} ~ 290+ at T=1200-2500), driving
# sigma_t to numerical zero and H_hat to garbage/NaN -- exactly the same
# failure mode diagnosed and fixed for the "fBM_persist" scenario in the
# ESN hyperparameter-search protocol. TIME_DT converts day-count to years.
TIME_DT = 1.0 / TRADING_DAYS

_CHOL_CACHE={}
def _get_chol(T,H,dt=1.0):
    """dt here is the ANNUALISATION factor (years/day-count-unit),
    default TIME_DT=1/252; NOT the simulation step size (always 1 day
    in dgp_lnsv's own loop)."""
    key=(T,H,dt)
    if key not in _CHOL_CACHE:
        t=np.arange(1,T+1,dtype=float)*dt
        C=0.5*(t[:,None]**(2*H)+t[None,:]**(2*H)-np.abs(t[:,None]-t[None,:])**(2*H))
        C+=np.eye(T)*1e-10
        try: L=sp_chol(C,lower=True)
        except: L=np.linalg.cholesky(C+np.eye(T)*1e-8)
        _CHOL_CACHE[key]=L
    return _CHOL_CACHE[key]

def dgp_lnsv(H,lam,n_paths,T,burn,sigma0=None,dt=1.0,seed_base=0,max_lag=40,q_grid=None):
    """
    r_t = sigma_t * eps_t,  eps_t ~ N(0,1)
    sigma_t = sigma_0 * exp( lambda * B^H(t) - lambda^2/2 * t^{2H} )
    Centring ensures E[sigma_t]=sigma_0.
    Cholesky simulation of B^H, with t measured in YEARS (TIME_DT),
    cached per (T,H,TIME_DT). `dt` here is the daily simulation step
    (always 1 trading day; kept as an argument for interface
    compatibility, not used to scale the fBM time argument).
    """
    if sigma0 is None: sigma0=TV_DAY
    if q_grid is None: q_grid=Q_GRID
    L=_get_chol(T,H,TIME_DT)
    rng_base=np.random.default_rng(seed_base); st_all=[]
    for p in range(n_paths):
        rng=np.random.default_rng(int(rng_base.integers(1<<31)))
        BH=L@rng.standard_normal(T)
        t_=np.arange(1,T+1,dtype=float)*TIME_DT
        log_sig=math.log(sigma0)+lam*BH-0.5*lam**2*t_**(2*H)
        sig=np.exp(log_sig); v=sig**2
        x   = sig * rng.standard_normal(T)          # full-length returns
        xb  = x[burn:]; vb = v[burn:]
        lv  = np.log(np.maximum(vb, 1e-30))
        # roughness diagnostics (for figures)
        st  = compute_roughness_stats(lv, max_lag, q_grid)
        # full 11-stat dict (for global score calibration)
        st_full = compute_statistics(xb, vb)
        st.update(st_full)                           # merge both stat dicts
        st["logret"] = xb; st["logV"] = lv
        st_all.append(st)
    return st_all

# ── Models (UNCHANGED) ────────────────────────────────────────
def _kernel_nodes(n_r,H,lam_lo=1/3500,lam_hi=2.0):
    lam_lo=max(float(lam_lo),1e-8); lam_hi=max(float(lam_hi),lam_lo*1.5)
    lam=np.geomspace(lam_lo,lam_hi,int(n_r))
    q=lam**(0.5-H); q/=(np.linalg.norm(q)+1e-15)
    b=np.sqrt(2.*lam)
    C=(b[:,None]*b[None,:])/(lam[:,None]+lam[None,:])
    q/=math.sqrt(float(q@C@q)+1e-15)
    return lam,q

def _build_esn(arch, H=None, lam_lo=1/3500, lam_hi=2.0, az_lo=1/280, az_hi=1/7,
               rough_scale=0.40, zr_lo=0.03, zr_hi=0.07, m1=0.0, m2=0.0):
    """
    Build reservoir matrices. H, lam_lo, lam_hi, az_lo, az_hi, rough_scale,
    zr_lo, zr_hi, m1, m2 are the per-(H,lambda)-cell CALIBRATED
    inner-loop parameters; if H is omitted, falls back to
    arch['H_target'] (used only for the start-up admissibility check).
    zz/sz/fi/si depend only on `arch` (Nr,Nz,zz_scale,sign_prob_neg,
    matrix_seed), not on H/lam/az/rough_scale/zr_lo/zr_hi/m1/m2.

    The read-out generalises z_readout to a per-mode profile
    z_readout_j = geomspace(zr_lo, zr_hi, Nz)[j], and adds a quadratic
    term r_t^T Q r_t, Q=(1/Nr)(m1*I+m2*qq^T) -- a genuine SYMMETRIC
    matrix (no positive-definiteness constraint) -- to the
    pre-activation (see _sim_esn). rough_scale, the z-readout profile,
    and (m1,m2) are all calibrated per scenario, not fixed architecture.
    """
    if H is None: H = arch["H_target"]
    rng=np.random.default_rng(int(arch["matrix_seed"]))
    n_r,n_z=arch["n_r"],arch["n_z"]
    lam,q=_kernel_nodes(n_r,H,lam_lo,lam_hi)
    az_lo=max(float(az_lo),1e-8); az_hi=max(float(az_hi),az_lo*1.5)
    az=np.geomspace(az_lo,az_hi,n_z)
    zr_lo_c = max(float(zr_lo), 1e-4); zr_hi_c = max(float(zr_hi), zr_lo_c*1.05)
    z_readout_profile = np.geomspace(zr_lo_c, zr_hi_c, n_z)
    zz=rng.uniform(-arch["zz_scale"],arch["zz_scale"],n_z)
    sz=-rng.choice([-1.,1.],n_z,p=[arch["sign_prob_neg"],1-arch["sign_prob_neg"]])
    fi=np.array([min(n_r-1,n_r//2+int((n_r//2-1)*j/max(n_z-1,1))) for j in range(n_z)],dtype=int)
    si=np.array([int((n_r//2-1)*(n_z-1-j)/max(n_z-1,1)) for j in range(n_z)],dtype=int)
    b0=_inv_sp(math.sqrt(max(TV_DAY**2-SIG_MIN**2,1e-15)))
    k0=arch["rough_orientation"]*rough_scale*float(q@np.sqrt(2.*lam))
    return dict(lam=lam,q=q,az=az,zz=zz,sz=sz,fi=fi,si=si,b0=b0,k0=k0,
                rough_scale=float(rough_scale),
                z_readout_profile=z_readout_profile,
                m1=float(m1), m2=float(m2))

_ESN_PARAMS=None

def _sim_esn(seed,T,dt=1.0,arch=_ARCH,params=None):
    """
    eta_t = b0 + rho_r*rough_scale*(q^T r_t) + sum_j w_{j,z} z_{j,t}
          + r_t^T Q r_t,
    Q = (1/Nr)*(m1*I + m2*qq^T), w_{j,z} = z_readout_j/sqrt(Nz),
    z_readout_j = geomspace(zr_lo, zr_hi, Nz)[j].
    rough_scale, the z_readout_profile, and (m1,m2) are read from
    `params` (calibrated per DGP cell), not from `arch`.
    """
    global _ESN_PARAMS
    if params is None:
        if _ESN_PARAMS is None: _ESN_PARAMS=_build_esn(arch)
        params=_ESN_PARAMS
    P=params; rng=np.random.default_rng(int(seed))
    n_r,n_z=arch["n_r"],arch["n_z"]
    spd=int(round(1./dt)); n_st=T*spd; sdt=math.sqrt(dt)
    al=np.exp(-P["lam"]*dt); cl=np.sqrt(np.maximum(1.-al**2,1e-14))
    azd=np.exp(-P["az"]*dt); om=1.-azd
    r=np.zeros(n_r); z=np.zeros(n_z)
    dx=np.zeros(T); dv=np.zeros(T)
    wz_vec = P["z_readout_profile"]/math.sqrt(n_z)
    rc=arch["rough_orientation"]*P["rough_scale"]
    m1c=P["m1"]/n_r; m2c=P["m2"]/n_r
    qvec=P["q"]
    for step in range(n_st):
        eps=rng.normal()
        qr=float(qvec@r)
        quad_term = m1c*float(r@r) + m2c*qr*qr
        eta=P["b0"]+rc*qr+float(wz_vec@z)+quad_term
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
            u=(P["sz"][j]*p1*p2+arch["even_strength"]*ev-arch["linear_strength"]*lin+P["zz"][j]*zo[j]+lc)
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

def sim_model_paths(model_name,cal,n_paths,T,burn,dt=1.0,max_lag=40,q_grid=None):
    if q_grid is None: q_grid=Q_GRID
    st_all=[]
    esn_params = None
    if model_name=="ESN A2-981003":
        esn_params = _build_esn(_ARCH, H=cal["H"], lam_lo=cal["lam_lo"], lam_hi=cal["lam_hi"],
                                 az_lo=cal["az_lo"], az_hi=cal["az_hi"],
                                 rough_scale=cal.get("rough_scale",0.40),
                                 zr_lo=cal.get("zr_lo",0.03), zr_hi=cal.get("zr_hi",0.07),
                                 m1=cal.get("m1",0.0), m2=cal.get("m2",0.0))
        esn_params["b0"] = esn_params["b0"] + cal.get("b0_delta",0.)
    for p in range(n_paths):
        if model_name=="ESN A2-981003":
            x,v=_sim_esn(9000+p,T,dt,arch=_ARCH,params=esn_params); x=x*cal.get("scale",1.0)
        elif model_name=="QRH":
            x,v=_sim_qrh(7000+p,T,dt,cal["H"],cal["nu_vol"],cal["lam"],cal["c_frac"])
        elif model_name=="PDV-GL":
            x,v=_sim_gl(8000+p,T,dt,cal["beta"],cal["alpha1"],cal["alpha_r"],cal["alpha_v"],cal["Vbar"])
        xb = x[burn:]; vb = v[burn:]
        lv = np.log(np.maximum(vb, 1e-30))
        st = compute_roughness_stats(lv, max_lag, q_grid)   # roughness diagnostics
        st.update(compute_statistics(xb, vb))               # full 11 score stats
        st["logret"] = xb; st["logV"] = lv
        st_all.append(st)
    return st_all

# ── Calibration (UNCHANGED) ───────────────────────────────────
def quick_calibrate(dgp_sts, n_cal=4, T_cal=600, burn_cal=150, dt=1.0):
    """
    Calibrate all three models using the FULL 11-term data-adaptive score S
    (smooth Gaussian-kernel version, minimised directly; S=0 is perfect).
    ref is built from the DGP cross-path statistics via make_score_ref().
    ESN calibrates its full 12-dim inner-loop vector (H, lam_lo, lam_hi,
    az_lo, az_hi, b0_delta, scale, rough_scale, zr_lo, zr_hi, m1, m2)
    via 11-start Nelder-Mead; only the 7 shared architecture
    hyperparameters (_ARCH) are held fixed.
    """

    # Build data-adaptive reference from all 11 DGP statistics
    ref = make_score_ref(dgp_sts)
    tH  = ref["H_hat_c"]
    print(f"    Global score ref: H={tH:.4f}  vol={ref['mean_vol_ann_c']*100:.1f}%  "
          f"Z={ref['zumbach_c']:.4f}  K={ref['kurtosis_c']:.3f}")

    # ── ESN: full inner-loop calibration (H, lam_lo, lam_hi, az_lo, az_hi, ─────
    #    b0_delta, scale, rough_scale, zr_lo, zr_hi, m1, m2), 11-point
    #    multi-start Nelder-Mead. Only the 7 shared hyperparameters of
    #    _ARCH are held fixed; rough_scale, the z-readout profile
    #    (zr_lo,zr_hi), and the symmetric quadratic matrix coefficients
    #    (m1,m2) are ALL calibrated per DGP cell here.
    rs_lo, rs_hi = ROUGH_SCALE_BOUNDS
    zrl_lo, zrl_hi = ZR_LO_BOUNDS
    zrh_lo, zrh_hi = ZR_HI_BOUNDS
    m1_lo, m1_hi = M1_BOUNDS
    m2_lo, m2_hi = M2_BOUNDS

    def esn_obj(p):
        H      = float(np.clip(p[0], 0.02, 0.49))
        lam_lo = float(np.exp(np.clip(p[1], -11, 1)))
        lam_hi = float(np.exp(np.clip(p[2], -6.0, 2.3)))
        az_lo  = float(np.exp(np.clip(p[3], -9, 0)))
        az_hi  = float(np.exp(np.clip(p[4], -5, 0)))
        b0d    = float(p[5])
        sc     = max(float(np.exp(np.clip(p[6], -1, 1))), 0.05)
        rough_scale = float(np.clip(math.exp(p[7]), rs_lo, rs_hi))
        zr_lo  = float(np.clip(math.exp(p[8]), zrl_lo, zrl_hi))
        zr_hi  = float(np.clip(math.exp(p[9]), zrh_lo, zrh_hi))
        m1     = float(np.clip(p[10], m1_lo, m1_hi))
        m2     = float(np.clip(p[11], m2_lo, m2_hi))
        lam_hi = max(lam_hi, lam_lo*1.5); az_hi = max(az_hi, az_lo*1.5)
        try:
            params = _build_esn(_ARCH, H=H, lam_lo=lam_lo, lam_hi=lam_hi, az_lo=az_lo, az_hi=az_hi,
                                 rough_scale=rough_scale, zr_lo=zr_lo, zr_hi=zr_hi,
                                 m1=m1, m2=m2)
            if params["k0"] >= 0: return SCORE_MAX
            params["b0"] = params["b0"] + b0d
            scores = []
            for q in range(n_cal):
                x, v = _sim_esn(3000+q, T_cal, dt, arch=_ARCH, params=params)
                xb = (x * sc)[burn_cal:]; vb = v[burn_cal:]
                scores.append(score_fn(compute_statistics(xb, vb), ref))
            return float(np.nanmean(scores))
        except Exception:
            return SCORE_MAX

    ESN_INNER_STARTS = [
        # (H0, lam_lo0, lam_hi0, az_lo0, az_hi0, rough_scale0, zr_lo0, zr_hi0, m1_0, m2_0)
        (0.10, 1/3500, 2.0,   1/280, 1/7,  0.40, 0.03, 0.07, 0.0, 0.0),  # rough / fast reservoir
        (0.40, 1e-4,   0.08,  1/280, 1/7,  0.40, 0.03, 0.07, 0.0, 0.0),  # persistent / slow reservoir
        (0.25, 1e-3,   0.5,   1/280, 1/7,  0.40, 0.03, 0.07, 0.0, 0.0),  # intermediate
        (0.20, 1e-5,   5.0,   1/400, 1/7,  0.40, 0.03, 0.07, 0.0, 0.0),  # wide-band / multiscale
        (0.35, 1e-5,   0.02,  1/500, 1/10, 0.40, 0.03, 0.07, 0.0, 0.0),  # long-memory extreme
        (0.15, 1e-4,   1.0,   1/300, 1/8,  0.90, 0.10, 0.30, 0.0, 0.0),  # strong amplitude
        (0.20, 1e-4,   1.0,   1/300, 1/8,  0.15, 0.02, 0.03, 0.0, 0.0),  # weak amplitude
        (0.20, 1e-4,   1.0,   1/300, 1/8,  0.40, 0.02, 0.20, 0.8, 0.0),  # isotropic +
        (0.20, 1e-4,   1.0,   1/300, 1/8,  0.40, 0.20, 0.02, -0.8, 0.0), # isotropic -
        (0.20, 1e-4,   1.0,   1/300, 1/8,  0.40, 0.03, 0.07, 0.0, 0.8),  # q-aligned +
        (0.20, 1e-4,   1.0,   1/300, 1/8,  0.40, 0.03, 0.07, 0.0, -0.8), # q-aligned -
    ]
    best_fun, best_x = None, None
    n_starts = len(ESN_INNER_STARTS)
    maxiter_per_start = max(30, 330 // n_starts)
    for (H0, ll0, lh0, al0, ah0, rs0, zrl0, zrh0, m10, m20) in ESN_INNER_STARTS:
        x0 = [H0, math.log(ll0), math.log(lh0), math.log(al0), math.log(ah0), 0.0, 0.0,
              math.log(rs0), math.log(zrl0), math.log(zrh0), m10, m20]
        res = optimize.minimize(esn_obj, x0, method="Nelder-Mead",
                                options={"maxiter": maxiter_per_start, "xatol": 0.02, "fatol": 0.02})
        if best_fun is None or res.fun < best_fun:
            best_fun, best_x = res.fun, res.x
    H      = float(np.clip(best_x[0], 0.02, 0.49))
    lam_lo = float(np.exp(np.clip(best_x[1], -11, 1)))
    lam_hi = float(np.exp(np.clip(best_x[2], -6.0, 2.3)))
    az_lo  = float(np.exp(np.clip(best_x[3], -9, 0)))
    az_hi  = float(np.exp(np.clip(best_x[4], -5, 0)))
    b0d    = float(best_x[5])
    sc     = max(float(np.exp(np.clip(best_x[6], -1, 1))), 0.05)
    rough_scale = float(np.clip(math.exp(best_x[7]), rs_lo, rs_hi))
    zr_lo  = float(np.clip(math.exp(best_x[8]), zrl_lo, zrl_hi))
    zr_hi  = float(np.clip(math.exp(best_x[9]), zrh_lo, zrh_hi))
    m1     = float(np.clip(best_x[10], m1_lo, m1_hi))
    m2     = float(np.clip(best_x[11], m2_lo, m2_hi))
    lam_hi = max(lam_hi, lam_lo*1.5); az_hi = max(az_hi, az_lo*1.5)
    print(f"    ESN: H={H:.4f}  lam=[{lam_lo:.2e},{lam_hi:.3f}]  "
          f"a=[{az_lo:.2e},{az_hi:.3f}]  b0_delta={b0d:.4f}  scale={sc:.4f}  "
          f"rough_scale={rough_scale:.4f}  zr=[{zr_lo:.4f},{zr_hi:.4f}]  "
          f"m1={m1:.4f}  m2={m2:.4f}  score={best_fun:.3f}")
    esn_cal = dict(H=H, lam_lo=lam_lo, lam_hi=lam_hi, az_lo=az_lo, az_hi=az_hi,
                   b0_delta=b0d, scale=sc, rough_scale=rough_scale,
                   zr_lo=zr_lo, zr_hi=zr_hi, m1=m1, m2=m2)

    # ── QRH ──────────────────────────────────────────────────────────────────
    def qrh_obj(p):
        H, nu, lk, cf = p
        if not (0.01<H<0.49 and 0.01<nu<3. and 0.1<lk<20. and 0.01<cf<0.99): return SCORE_MAX
        scores = []
        for q in range(n_cal):
            try:
                x, v = _sim_qrh(500+q, T_cal, dt, H, nu, lk, cf)
                scores.append(score_fn(compute_statistics(x[burn_cal:], v[burn_cal:]), ref))
            except: scores.append(SCORE_MAX)
        return float(np.nanmean(scores))
    res = optimize.minimize(qrh_obj,
                            [max(0.02, min(tH, 0.48)), 0.30, 2.0, 0.50],
                            method="Nelder-Mead",
                            options={"maxiter": 300, "xatol": 0.01, "fatol": 0.02})
    H_o, nv_o, lk_o, cf_o = res.x
    H_o  = float(np.clip(H_o,  0.01, 0.48)); nv_o = float(np.clip(nv_o, 0.01, 3.))
    lk_o = float(np.clip(lk_o, 0.1,  20.)); cf_o = float(np.clip(cf_o, 0.01, 0.99))
    print(f"    QRH: H={H_o:.4f}  nu={nv_o:.4f}  lam={lk_o:.4f}  c={cf_o:.4f}  score={res.fun:.3f}")
    qrh_cal = dict(H=H_o, nu_vol=nv_o, lam=lk_o, c_frac=cf_o)

    # ── PDV-GL ───────────────────────────────────────────────────────────────
    Vbar = TV_DAY**2 * dt
    def gl_obj(p):
        beta, a1, ar, av = p
        if not (0.05<beta<0.99 and 0<a1<1 and 0.01<ar<0.49 and 0.01<av<0.49): return SCORE_MAX
        scores = []
        for q in range(n_cal):
            try:
                x, v = _sim_gl(600+q, T_cal, dt, beta, a1, ar, av, Vbar)
                scores.append(score_fn(compute_statistics(x[burn_cal:], v[burn_cal:]), ref))
            except: scores.append(SCORE_MAX)
        return float(np.nanmean(scores))
    res = optimize.minimize(gl_obj, [0.75, 0.45, 0.08, 0.25], method="Nelder-Mead",
                            options={"maxiter": 400, "xatol": 0.01, "fatol": 0.02})
    b_o, a1_o, ar_o, av_o = res.x
    b_o  = float(np.clip(b_o,  0.05, 0.99)); a1_o = float(np.clip(a1_o, 0.01, 0.99))
    ar_o = float(np.clip(ar_o, 0.01, 0.49)); av_o = float(np.clip(av_o, 0.01, 0.49))
    print(f"    PDV-GL: beta={b_o:.4f}  a1={a1_o:.4f}  ar={ar_o:.4f}  av={av_o:.4f}  score={res.fun:.3f}")
    gl_cal = dict(beta=b_o, alpha1=a1_o, alpha_r=ar_o, alpha_v=av_o, Vbar=Vbar)

    return esn_cal, qrh_cal, gl_cal, ref   # return ref so run() can store it

# ── Figures ───────────────────────────────────────────────────
def _violin(ax,data,pos,color,width=0.35):
    data=np.array([d for d in data if np.isfinite(d)])
    if len(data)<3: return
    try:
        parts=ax.violinplot([data],positions=[pos],widths=width,
                             showmeans=False,showmedians=True,showextrema=False)
        for pc in parts["bodies"]: pc.set_facecolor(color); pc.set_alpha(0.45)
        parts["cmedians"].set_color(color); parts["cmedians"].set_linewidth(1.5)
    except: pass
    jit=np.random.default_rng(42).uniform(-width*0.3,width*0.3,len(data))
    ax.scatter(pos+jit,data,color=color,s=6,alpha=0.45,zorder=4)

def plot_hurst_estimation(results,H_grid,lam_fixed,save_path=None):
    Hcols=[H for H in H_grid if (H,lam_fixed) in results]
    n_c=len(Hcols); all_names=["DGP (fBM SV)"]+MODEL_NAMES
    fig,axes=plt.subplots(3,n_c,figsize=(5*n_c,13))
    if n_c==1: axes=axes[:,None]
    fig.suptitle(f"$\hat H$ Estimation — log-normal SV DGP ($\lambda={lam_fixed}$)\n"
                 r"Dyadic variogram OLS: slope $=2\hat H$",fontsize=11,fontweight="bold")
    for ci,H in enumerate(Hcols):
        data=results[(H,lam_fixed)]; axes[0,ci].set_title(f"True $H={H}$",fontsize=11,fontweight="bold")
        ax=axes[0,ci]
        for mi,name in enumerate(all_names):
            Hhats=[s["H_hat"] for s in data[name] if np.isfinite(s["H_hat"])]
            _violin(ax,Hhats,mi,COLORS.get(name,GRAY))
        ax.axhline(H,color=CORAL,lw=1.4,ls="--",zorder=6,label=f"True $H={H}$")
        dgp_Hs=[s["H_hat"] for s in data["DGP (fBM SV)"] if np.isfinite(s["H_hat"])]
        if len(dgp_Hs)>=10:
            lo=float(np.percentile(dgp_Hs,2.5)); hi=float(np.percentile(dgp_Hs,97.5))
            ax.axhspan(lo,hi,alpha=0.09,color=GRAY,label=f"DGP 95% CI")
        ax.set_xticks(range(len(all_names))); ax.set_xticklabels(all_names,rotation=25,fontsize=7,ha="right")
        if ci==0: ax.set_ylabel("$\hat H$",fontsize=10)
        ax.legend(fontsize=7); ax.grid(True,axis="y",alpha=0.18)
        ax=axes[1,ci]
        for mi,name in enumerate(all_names):
            Hhats=np.array([s["H_hat"] for s in data[name] if np.isfinite(s["H_hat"])])
            if len(Hhats)==0: continue
            bias=float(Hhats.mean())-H; se=float(Hhats.std())/math.sqrt(len(Hhats))
            ax.bar(mi,bias,color=COLORS.get(name,GRAY),alpha=0.82,width=0.7)
            ax.errorbar(mi,bias,yerr=1.96*se,color="k",capsize=4,lw=1.2)
        ax.axhline(0,color=GRAY,lw=0.9,ls="--")
        ax.set_xticks(range(len(all_names))); ax.set_xticklabels(all_names,rotation=25,fontsize=7,ha="right")
        if ci==0: ax.set_ylabel("Bias",fontsize=9)
        ax.set_title("Bias (±1.96 SE)",fontsize=9); ax.grid(True,axis="y",alpha=0.18)
        ax=axes[2,ci]
        for mi,name in enumerate(all_names):
            Hhats=np.array([s["H_hat"] for s in data[name] if np.isfinite(s["H_hat"])])
            if len(Hhats)==0: continue
            rmse=float(np.sqrt(np.mean((Hhats-H)**2)))
            ax.bar(mi,rmse,color=COLORS.get(name,GRAY),alpha=0.82,width=0.7)
            ax.text(mi,rmse+0.002,f"{rmse:.3f}",ha="center",fontsize=7)
        ax.set_xticks(range(len(all_names))); ax.set_xticklabels(all_names,rotation=25,fontsize=7,ha="right")
        if ci==0: ax.set_ylabel("RMSE$(\hat H)$",fontsize=9)
        ax.set_title("RMSE",fontsize=9); ax.grid(True,axis="y",alpha=0.18)
    plt.tight_layout(rect=[0,0,1,0.94])
    if save_path: fig.savefig(save_path,dpi=150,bbox_inches="tight"); print(f"  Saved -> {save_path}")
    return fig

def plot_acf_scaling(results,H_grid,lam_fixed,short_lag=20,save_path=None):
    Hcols=[H for H in H_grid if (H,lam_fixed) in results]
    all_names=["DGP (fBM SV)"]+MODEL_NAMES
    n_c=len(Hcols); n_r=len(all_names)
    fig,axes=plt.subplots(n_r,n_c,figsize=(5*n_c,4*n_r))
    if n_c==1: axes=axes[:,None]
    if n_r==1: axes=axes[None,:]
    fig.suptitle(f"ACov of $\log V_t$ vs $\tau^{{2\hat H}}$ — ($\lambda={lam_fixed}$)\n"
                 r"OLS: $C(\tau)/C(0)=a+b\tau^{2\hat H}$  |  $R^2\to1$ = correct roughness",
                 fontsize=11,fontweight="bold")
    for ci,H in enumerate(Hcols):
        data=results[(H,lam_fixed)]; axes[0,ci].set_title(f"True $H={H}$",fontsize=10,fontweight="bold")
        for ri,name in enumerate(all_names):
            ax=axes[ri,ci]; sts=data[name]
            C0s=np.array([s["C0"] for s in sts if s["C0"]>0])
            Cks=np.array([s["Ck"] for s in sts if s["C0"]>0])
            if len(C0s)==0: continue
            C0_m=float(C0s.mean()); Ck_m=Cks.mean(0)
            H_m=float(np.nanmean([s["H_hat"] for s in sts]))
            reg=acf_scaling_regression(C0_m,Ck_m,H_m,short_lag)
            ax.scatter(reg["tau_2H"],reg["rho"],color=COLORS.get(name,GRAY),s=7,alpha=0.55,zorder=3)
            if np.isfinite(reg["a"]) and np.isfinite(reg["b"]):
                xf=reg["tau_2H"][:reg["n_fit"]]
                ax.plot(xf,reg["a"]+reg["b"]*xf,color=CORAL,lw=1.8,ls="--",zorder=5,
                        label=f"OLS $R^2={reg['r2']:.3f}$")
            ax.axhline(0,color=GRAY,lw=0.5,ls="--",alpha=0.5)
            ax.annotate(f"$\hat H={H_m:.3f}$\n$a={reg['a']:.3f}$  $b={reg['b']:.3f}$\n$R^2={reg['r2']:.3f}$",
                        xy=(0.97,0.95),xycoords="axes fraction",ha="right",va="top",fontsize=7.5,
                        family="monospace",bbox=dict(boxstyle="round,pad=0.3",fc="white",ec="#ccc",alpha=0.92))
            if ci==0: ax.set_ylabel(f"{name}\n$C(\tau)/C(0)$",fontsize=8)
            if ri==n_r-1: ax.set_xlabel(r"$\tau^{2\hat H}$",fontsize=9)
            ax.legend(fontsize=7); ax.grid(True,alpha=0.18)
    plt.tight_layout(rect=[0,0,1,0.93])
    if save_path: fig.savefig(save_path,dpi=150,bbox_inches="tight"); print(f"  Saved -> {save_path}")
    return fig

def plot_scale_invariance(results,H_grid,lam_fixed,save_path=None):
    Hcols=[H for H in H_grid if (H,lam_fixed) in results]
    all_names=["DGP (fBM SV)"]+MODEL_NAMES
    n_c=len(all_names); q_arr=np.array(Q_GRID,float)
    marks=["o","s","^","D"]; lss=["-","--","-.",":"]
    fig,axes=plt.subplots(3,n_c,figsize=(5*n_c,13))
    if n_c==1: axes=axes[:,None]
    fig.suptitle(f"Scale Invariance — log-normal SV DGP ($\lambda={lam_fixed}$)\n"
                 r"$S_q(\tau)\sim\tau^{\zeta(q)}$  |  Monofractal DGP: $\zeta(q)=qH$",
                 fontsize=11,fontweight="bold")
    for ci,name in enumerate(all_names):
        col=COLORS.get(name,GRAY)
        for ri in range(3): axes[ri,ci].grid(True,alpha=0.18)
        axes[0,ci].set_title(name,fontsize=10,fontweight="bold")
        for hi,H in enumerate(Hcols):
            data=results[(H,lam_fixed)][name]
            zetas=np.array([s["zeta"] for s in data]); Hqs=np.array([s["Hq"] for s in data])
            mus=np.array([s["mu"] for s in data])
            zm=np.nanmean(zetas,axis=0); zs=np.nanstd(zetas,axis=0)
            Hm=np.nanmean(Hqs,axis=0); Hs=np.nanstd(Hqs,axis=0)
            mm=np.nanmean(mus,axis=0); ms=np.nanstd(mus,axis=0)
            lbl=f"$H={H}$"; mk=marks[hi%len(marks)]; ls=lss[hi%len(lss)]
            ax=axes[0,ci]
            ax.fill_between(q_arr,zm-zs,zm+zs,alpha=0.12,color=col)
            ax.plot(q_arr,zm,color=col,lw=1.8,marker=mk,ms=5,ls=ls,label=lbl)
            ax.plot(q_arr,q_arr*H,color=CORAL,lw=0.8,ls=":",alpha=0.6)
            ax=axes[1,ci]
            ax.fill_between(q_arr,Hm-Hs,Hm+Hs,alpha=0.12,color=col)
            ax.plot(q_arr,Hm,color=col,lw=1.8,marker=mk,ms=5,ls=ls,label=lbl)
            ax.axhline(H,color=CORAL,lw=0.8,ls=":",alpha=0.6)
            ax=axes[2,ci]
            ax.fill_between(q_arr,mm-ms,mm+ms,alpha=0.12,color=col)
            ax.plot(q_arr,mm,color=col,lw=1.8,marker=mk,ms=5,ls=ls,label=lbl)
        if ci==0:
            axes[0,ci].set_ylabel(r"$\zeta(q)$",fontsize=10)
            axes[1,ci].set_ylabel(r"$H(q)=\zeta(q)/q$",fontsize=10)
            axes[2,ci].set_ylabel(r"$\mu(q)=\zeta(q)-q\zeta(1)$",fontsize=10)
        for ri in range(3):
            axes[ri,ci].set_xlabel("$q$" if ri==2 else "",fontsize=9)
            axes[ri,ci].legend(fontsize=7)
        axes[2,ci].axhline(0,color=GRAY,lw=0.8,ls="--",alpha=0.6)
    plt.tight_layout(rect=[0,0,1,0.93])
    if save_path: fig.savefig(save_path,dpi=150,bbox_inches="tight"); print(f"  Saved -> {save_path}")
    return fig

def plot_hurst_summary(results,H_grid,lam_fixed,save_path=None):
    Hcols=[H for H in H_grid if (H,lam_fixed) in results]
    all_names=["DGP (fBM SV)"]+MODEL_NAMES
    all_cols=[COLORS.get(n,GRAY) for n in all_names]; marks=["o","s","^","D"]
    fig,axes=plt.subplots(3,1,figsize=(10,13))
    fig.suptitle(f"Roughness Summary vs True $H$ — ($\lambda={lam_fixed}$)\n"
                 r"$\hat H$ via OLS dyadic variogram  |  Bands = 95% cross-path CI",
                 fontsize=11,fontweight="bold")
    for mi,(name,col,mk) in enumerate(zip(all_names,all_cols,marks)):
        Hm_l,lo_l,hi_l,bi_l,rm_l=[],[],[],[],[]
        for H in Hcols:
            Hhats=np.array([s["H_hat"] for s in results[(H,lam_fixed)][name] if np.isfinite(s["H_hat"])])
            if len(Hhats)==0:
                for lst in [Hm_l,lo_l,hi_l,bi_l,rm_l]: lst.append(np.nan)
                continue
            Hm_l.append(float(Hhats.mean())); lo_l.append(float(np.percentile(Hhats,2.5)))
            hi_l.append(float(np.percentile(Hhats,97.5))); bi_l.append(float(Hhats.mean())-H)
            rm_l.append(float(np.sqrt(np.mean((Hhats-H)**2))))
        Hc=np.array(Hcols); Hm_a=np.array(Hm_l); lo_a=np.array(lo_l); hi_a=np.array(hi_l)
        off=(mi-len(all_names)/2)*0.007
        axes[0].errorbar(Hc+off,Hm_a,yerr=[Hm_a-lo_a,hi_a-Hm_a],
                         fmt=mk,color=col,ms=7,lw=1.5,capsize=4,label=name)
        axes[1].plot(Hc+off,bi_l,mk+"-",color=col,ms=7,lw=1.5,label=name)
        axes[2].plot(Hc+off,rm_l,mk+"-",color=col,ms=7,lw=1.5,label=name)
    Hr=np.linspace(min(Hcols)*0.8,max(Hcols)*1.1,100)
    axes[0].plot(Hr,Hr,color=CORAL,lw=1.2,ls="--",label="Identity $\hat H=H$")
    axes[1].axhline(0,color=CORAL,lw=1.2,ls="--",label="Zero bias")
    for ri,(ax,yl,tt) in enumerate(zip(axes,[
            r"$\bar{\hat H}$ with 95% CI",r"Bias $=\bar{\hat H}-H$",r"RMSE$(\hat H)$"],[
            r"$\hat H$ vs true $H$  (identity=perfect)",r"Bias (zero=unbiased)",r"RMSE"])):
        ax.set_xticks(Hcols); ax.set_xticklabels([f"${H}$" for H in Hcols])
        ax.set_xlabel("True $H$" if ri==2 else "",fontsize=10)
        ax.set_ylabel(yl,fontsize=10); ax.set_title(tt,fontsize=10)
        ax.legend(fontsize=8,ncol=2); ax.grid(True,alpha=0.20)
    plt.tight_layout(rect=[0,0,1,0.95])
    if save_path: fig.savefig(save_path,dpi=150,bbox_inches="tight"); print(f"  Saved -> {save_path}")
    return fig

# ── LaTeX ─────────────────────────────────────────────────────
LATEX = r"""\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb,amsthm,amsfonts}
\usepackage{geometry}\geometry{left=2.3cm,right=2.3cm,top=2.5cm,bottom=2.5cm}
\usepackage{setspace}\setstretch{1.20}
\usepackage{hyperref}\hypersetup{colorlinks=true,linkcolor=black,citecolor=black}
\usepackage{parskip,booktabs,array,xcolor,titlesec,enumitem,microtype}
\microtypesetup{expansion=false}
\titleformat{\section}{\large\bfseries}{}{0em}{}[\titlerule]
\titlespacing{\section}{0pt}{14pt}{6pt}
\titleformat{\subsection}{\normalsize\bfseries}{\thesubsection\;}{0em}{}
\newtheorem{remark}{Remark}
\begin{document}
\begin{center}
{\LARGE\bfseries Numerical Protocol}\\[0.3cm]
{\large Log-Volatility Roughness Evaluation\\
under the Log-Normal SV DGP with fBM}\\[0.2cm]
{\normalsize Othmane Zarhali --- Paris Dauphine / CNRS}
\end{center}
\vspace{0.15cm}\noindent\rule{\textwidth}{0.8pt}
\begin{abstract}
We evaluate the roughness of $\log V_t$ under the log-normal SV DGP
$\sigma_t=\sigma_0\exp(\lambda B^H(t)-\lambda^2 t^{2H}/2)$,
simulated via Cholesky decomposition of the fBM covariance matrix.
The DGP has analytically exact variogram slope $2H$ and is monofractal:
$\zeta(q)=qH$, $H(q)\equiv H$, $\mu(q)\equiv0$.
Three models (ESN A2-981003, run at its current optimal architecture,
$N_r{=}96,N_z{=}12$, using the \textbf{smooth Gaussian-kernel score} --
QRH, and PDV-GL) are calibrated to the DGP's full stylised-fact profile
and evaluated on three roughness diagnostics: OLS Hurst estimation,
ACov vs $\tau^{2\hat H}$ linearity, and multifractal spectrum.
\end{abstract}
\tableofcontents\vspace{0.4cm}
\section{DGP: Log-Normal SV with fBM}
\begin{equation}
  r_t=\sigma_t\varepsilon_t,\quad\varepsilon_t\sim\mathcal{N}(0,1),\qquad
  \sigma_t=\sigma_0\exp\!\Bigl(\lambda B^H(t)-\tfrac{\lambda^2}{2}t^{2H}\Bigr).
\end{equation}
The centring term $-\lambda^2 t^{2H}/2$ ensures $\mathbb{E}[\sigma_t]=\sigma_0$.
Analytical truth: $m_2(\tau)=4\lambda^2\tau^{2H}$, so OLS variogram slope $=2H$ exactly.
DGP is monofractal: $\zeta(q)=qH$, $H(q)\equiv H$, $\mu(q)\equiv0$.
Simulation via Cholesky of $\Sigma_{ij}=\tfrac12(t_i^{2H}+t_j^{2H}-|t_i-t_j|^{2H})$, cached.
Parameter grid available: $H\in\{0.05,0.10,0.40\}$, $\lambda\in\{0.5,1.0,2.0\}$; this
note reports the rough-target subset $H\in\{0.05,0.10\}$ at $\lambda=1.0$.
\section{Models}
All three models use Gaussian Brownian innovations.
\subsection{ESN A2-981003 (current optimal architecture)}
7 shared hyperparameters, with $N_r=96,N_z=12$.
\texttt{rough\_scale}, the $z$-readout profile ($\texttt{zr\_lo},
\texttt{zr\_hi}$), and the symmetric quadratic matrix coefficients
$(m_1,m_2)$ are calibrated per $(H,\lambda)$ cell, not part of this
fixed table (see \S3's Calibration table for their per-cell values).
\begin{table}[h]
\centering
\caption{The 7 shared hyperparameters of the optimal architecture.}
\begin{tabular}{@{}lr@{}}
\toprule
Hyperparameter & Value \\
\midrule
$N_r$ & 96 \\
$N_z$ & 12 \\
\texttt{z\_strength} & 0.1929 \\
\texttt{even\_strength} & 2.2291 \\
\texttt{linear\_strength} & 0.1930 \\
\texttt{gamma\_norm} & 0.9651 \\
\texttt{local\_z\_strength} & 0.0476 \\
\texttt{zz\_scale} & 0.0397 \\
\texttt{sign\_prob\_neg} & 0.2219 \\
\texttt{rough\_orientation} & $-1.0$ (fixed) \\
\bottomrule
\end{tabular}
\end{table}
A single $\varepsilon_t\sim\mathcal{N}(0,1)$ drives both the return and
the reservoir. $H$, the $r$-layer rate window
$[\lambda_{\rm lo},\lambda_{\rm hi}]$, and the $z$-layer rate window
$[a_{\rm lo},a_{\rm hi}]$ are \textbf{calibrated per $(H,\lambda)$ DGP
cell} (\S3) -- they are not part of the shared architecture above.
$H_{\rm ESN}$ (the reservoir's own memory-kernel exponent) is a
different object from the DGP's Hurst parameter $H_{\rm DGP}$; the
inner loop searches over $H_{\rm ESN}$ so the simulated $\hat H$
matches the DGP's target, and the two need not coincide.
\subsection{QRH}
Bourgey \& Gatheral (2026): $V_t=Y_t^2+c$, gamma kernel, Euler-Volterra
with RL weights $w_j=\Delta t^\alpha/(\Gamma(\alpha)\alpha)[(j+1)^\alpha-j^\alpha]$,
252-step ring buffer.
\subsection{PDV-GL}
Guyon \& Lekeufack (2023): two-factor power-law kernel,
$V_t=(1-\beta)\bar V+\beta[\alpha_1 F_t(\alpha_r)+(1-\alpha_1)F_t(\alpha_v)]$,
online computation.
\section{Calibration}
Calibration minimises the full 11-term data-adaptive score $S$, with
the \textbf{smooth Gaussian-kernel} proximity:
\begin{equation}
S=\sum_k w_k\bigl(1-f_k\bigr)+5\cdot\mathrm{Stress},\qquad
f(x,c,s)=\exp\!\Bigl(-\tfrac12\bigl(\tfrac{x-c}{s}\bigr)^2\Bigr)\in(0,1],
\end{equation}
$S\geq0$, $S=0$ is a perfect match on all 11 stylised facts, $S=15.3$ is
worst-possible; minimised directly (no sign flip). The smooth kernel
provides a never-zero gradient at any distance from target, unlike a
piecewise-linear tent $\max(0,1-|x-c|/s)$, which hits exactly 0 (zero
gradient) beyond one tolerance width.

\paragraph{ESN read-out.} The ESN's volatility pre-activation combines
three calibrated channels:
\begin{equation}
\eta_t = b_0 + \rho_r\,\texttt{rough\_scale}\,(q^\top r_t)
       + \sum_{j=1}^{N_z} w_{j,z}\,z_{j,t} + r_t^\top Qr_t,
\qquad w_{j,z}=\frac{\texttt{z\_readout}_j}{\sqrt{N_z}},
\qquad Q=\frac{1}{N_r}(m_1I+m_2qq^\top),
\end{equation}
where $\texttt{z\_readout}_j=\mathrm{geomspace}(\texttt{zr\_lo},
\texttt{zr\_hi},N_z)[j]$ is a per-mode profile, and $Q$ is a genuine
\textbf{symmetric} matrix (no positive-definiteness constraint): $m_1$
scales the reservoir's isotropic energy, $m_2$ scales the squared rough
factor $(q^\top r_t)^2$ specifically. \textbf{Only 7 hyperparameters
remain shared} across all $(H,\lambda)$ cells
(\texttt{z\_strength}, \texttt{even\_strength}, \texttt{linear\_strength},
\texttt{gamma\_norm}, \texttt{local\_z\_strength}, \texttt{zz\_scale},
\texttt{sign\_prob\_neg}, plus $N_r,N_z$); \texttt{rough\_scale},
\texttt{zr\_lo}, \texttt{zr\_hi}, and $(m_1,m_2)$ are calibrated
per $(H,\lambda)$ cell, exactly as $H,\lambda_{\rm lo},\lambda_{\rm hi},
a_{\rm lo},a_{\rm hi}$ already were.

\begin{center}\renewcommand{\arraystretch}{1.3}
\begin{tabular}{@{}p{3.0cm}p{6.4cm}p{5.0cm}@{}}\toprule
Model & Parameters calibrated & Remark \\\midrule
ESN A2-981003 & $H$, $\lambda_{\rm lo}$, $\lambda_{\rm hi}$, $a_{\rm lo}$,
  $a_{\rm hi}$, $\Delta b_0$, scale $s$, \texttt{rough\_scale},
  \texttt{zr\_lo}, \texttt{zr\_hi}, $m_1$, $m_2$ &
  7 shared hyperparameters fixed; full 12-dim.\ inner-loop calibrated
  per $(H,\lambda)$ cell (11-start Nelder--Mead) \\
QRH & $H$, $\nu_{\rm vol}$, $\lambda$, $c_{\rm frac}$ & \\
PDV-GL & $\beta$, $\alpha_1$, $\alpha_r$, $\alpha_v$ & \\
\bottomrule\end{tabular}\end{center}
The ESN's inner loop uses an 11-point multi-start Nelder--Mead scheme
(rough/fast, persistent/slow, intermediate, wide-band/multiscale,
long-memory extreme, strong-amplitude, weak-amplitude,
isotropic-positive, isotropic-negative, $q$-aligned-positive, and
$q$-aligned-negative starts), keeping the best-of-eleven result.
\section{Per-Cell Calibrated Read-Out Parameters}
At $\lambda=1.0$, the ESN's calibrated \texttt{rough\_scale} and
$z$-readout profile stay in a narrow band across the two rough-target
$H$ values tested; $m_1,m_2$ are both near-zero at both cells, with no
large excursion in either quadratic coefficient at this narrower $H$
range (this contrasts with a persistent target, $H=0.40$, tested
separately, where $m_1$ calibrated to a large negative value --
consistent with the general caution that this quadratic channel is
weakly and inconsistently identified across seeds and should not be
over-interpreted from any single cell).
\begin{center}\renewcommand{\arraystretch}{1.3}
\begin{tabular}{@{}lcccccc@{}}\toprule
True $H$ & \texttt{rough\_scale} & \texttt{zr\_lo} & \texttt{zr\_hi} & $m_1$ & $m_2$ & Score $S$ \\\midrule
0.05 & 0.433 & 0.032 & 0.072 & $-0.0002$ & $-0.0002$ & 3.955 \\
0.10 & 0.410 & 0.028 & 0.061 & $0.0000$ & $0.0004$ & 4.216 \\
\bottomrule\end{tabular}\end{center}

\section{Roughness Diagnostics}
\subsection{Diagnostic 1: Hurst OLS}
$\hat H=\tfrac12\cdot\operatorname{OLS slope}$ of $\log m_2(L)$ vs $\log L$,
$L\in\{1,2,4,8,16,32,64\}$. Cross-path 95\% CI from empirical quantiles.
\subsection{Diagnostic 2: ACov linearity}
OLS: $C(\tau)/C(0)=a+b\tau^{2\hat H}$. $R^2\to1$ = correct Hölder exponent.
\subsection{Diagnostic 3: Scale invariance}
$S_q(\tau)\sim\tau^{\zeta(q)}$; $H(q)=\zeta(q)/q$; $\mu(q)=\zeta(q)-q\zeta(1)$.
Monofractal DGP: all three are constant/zero.
\section{Figures}
Per $\lambda$: R1 violin+bias+RMSE; R2 ACov scatter+OLS; R3 $\zeta,H(q),\mu(q)$; R4 summary.
\section{Monte Carlo Setup}
\begin{center}\renewcommand{\arraystretch}{1.4}
\begin{tabular}{ll}\toprule Parameter & Value \\\midrule
DGP paths & 30 \\ Cal.\ paths & 4 \\ Cal.\ days & 600 \\
Final paths & 30 \\ Final days & 2500 \\ Burn-in & 500 \\
$\Delta t$ & 1 (daily) \\ Dyadic lags & $\{1,2,4,8,16,32,64\}$ \\
Structure function $q$ & $\{0.5,1.0,1.5,2.0,2.5,3.0\}$ \\
ACov lags $\tau_{\max}$ & 20 \\
\bottomrule\end{tabular}\end{center}
\begin{thebibliography}{9}
\bibitem{gatheral2018} J.\ Gatheral, T.\ Jaisson, M.\ Rosenbaum.
\textit{Volatility is rough}. QF 18(6), 2018.
\bibitem{bourgey2026} F.\ Bourgey, J.\ Gatheral.
\textit{Quadratic Rough Heston}. SSRN:5239929, 2026.
\bibitem{guyon2023} J.\ Guyon, J.\ Lekeufack.
\textit{Volatility is (mostly) path-dependent}. QF 23(9), 2023.
\end{thebibliography}
\end{document}
"""

def write_latex(path):
    with open(path,"w") as f: f.write(LATEX.lstrip())
    print(f"  Saved -> {path}")

# ── Main pipeline ─────────────────────────────────────────────
def run(H_grid=None,lam_grid=None,lam_fixed=None,n_dgp=30,n_sim=30,
        T=2500,burn=500,n_cal=4,T_cal=600,burn_cal=150,
        dt=1.0,max_lag=40,q_grid=None,short_lag=20,
        save_prefix="roughness_lnsv",latex_path=None):
    if H_grid is None: H_grid=H_GRID
    if lam_grid is None: lam_grid=LAMBDA_GRID
    if lam_fixed is None: lam_fixed=lam_grid[min(1,len(lam_grid)-1)]
    if q_grid is None: q_grid=Q_GRID
    global _ESN_PARAMS; _ESN_PARAMS=_build_esn(_ARCH)
    assert _ESN_PARAMS["k0"]<0; print(f"ESN kappa0={_ESN_PARAMS['k0']:.6f}  ✓\n")
    results={}
    for lam in lam_grid:
        for H in H_grid:
            key=(H,lam); print(f"{'='*60}  H={H}  lambda={lam}")
            print("  DGP (log-normal SV fBM) ...")
            dgp_sts=dgp_lnsv(H,lam,n_dgp,T,burn,dt=dt,
                              seed_base=int(H*1000)+int(lam*100),
                              max_lag=max_lag,q_grid=q_grid)
            Hm=float(np.nanmean([s["H_hat"] for s in dgp_sts]))
            print(f"  DGP: H_hat={Hm:.4f}  (true H={H})")
            print("  Calibrating (global 11-term score S) ...")
            esn_cal, qrh_cal, gl_cal, ref = quick_calibrate(dgp_sts, n_cal, T_cal, burn_cal, dt)
            print("  Simulating final paths ...")
            esn_sts=sim_model_paths("ESN A2-981003",esn_cal,n_sim,T,burn,dt,max_lag,q_grid)
            qrh_sts=sim_model_paths("QRH",qrh_cal,n_sim,T,burn,dt,max_lag,q_grid)
            gl_sts=sim_model_paths("PDV-GL",gl_cal,n_sim,T,burn,dt,max_lag,q_grid)
            results[key] = {"DGP (fBM SV)": dgp_sts, "ESN A2-981003": esn_sts,
                            "QRH": qrh_sts, "PDV-GL": gl_sts, "ref": ref}
            print(f"\n  Results (true H={H}, lambda={lam}):")
            print(f"  {'Model':<22} {'H_hat':>7} {'vol%':>6} {'Z':>8} {'Lev':>7} {'Kurt':>7} {'Score':>7}")
            print("  " + "-"*62)
            for nm in ["DGP (fBM SV)"] + MODEL_NAMES:
                sts = results[key][nm]
                def _m(k): return float(np.nanmean([s[k] for s in sts if np.isfinite(s.get(k, np.nan))]))
                sc = float(np.nanmean([score_fn(s, ref) for s in sts]))
                print(f"  {nm:<22} {_m('H_hat'):>7.4f} "
                      f"{_m('mean_vol_ann')*100:>5.1f}% "
                      f"{_m('zumbach'):>8.4f} "
                      f"{_m('leverage'):>7.4f} "
                      f"{_m('kurtosis'):>7.3f} {sc:>7.3f}")
            print()
    print(f"\n{'='*60}\nFigures ...\n{'='*60}")
    for lam in lam_grid:
        H_ok=[H for H in H_grid if (H,lam) in results]
        if not H_ok: continue
        sfx=f"lam{str(lam).replace('.','p')}"
        for fn,sp in [(plot_hurst_estimation,"hurst"),(plot_acf_scaling,"acf_scaling"),
                      (plot_scale_invariance,"scale_inv"),(plot_hurst_summary,"hurst_summary")]:
            if fn==plot_acf_scaling:
                fig=fn(results,H_ok,lam,short_lag=short_lag,save_path=f"{save_prefix}_{sp}_{sfx}.png")
            else:
                fig=fn(results,H_ok,lam,save_path=f"{save_prefix}_{sp}_{sfx}.png")
            plt.close(fig)
    lp=latex_path or f"{save_prefix}_protocol.tex"
    write_latex(lp)
    if os.system("which pdflatex>/dev/null 2>&1")==0:
        os.system(f"pdflatex -interaction=nonstopmode {lp}>/dev/null 2>&1")
        os.system(f"pdflatex -interaction=nonstopmode {lp}>/dev/null 2>&1")
        pdf=lp.replace(".tex",".pdf")
        if os.path.exists(pdf): print(f"  Compiled -> {pdf}")
    return results

if __name__=="__main__":
    import argparse
    pa=argparse.ArgumentParser(description="Log-vol roughness under log-normal SV DGP",
                               formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pa.add_argument("--H",nargs="+",type=float,default=[0.05,0.10,0.40])
    pa.add_argument("--lam",nargs="+",type=float,default=[0.5,1.0,2.0])
    pa.add_argument("--lam_fixed",type=float,default=1.0)
    pa.add_argument("--n_dgp",type=int,default=30); pa.add_argument("--n_cal",type=int,default=4)
    pa.add_argument("--n_sim",type=int,default=30); pa.add_argument("--T",type=int,default=2500)
    pa.add_argument("--burn",type=int,default=500); pa.add_argument("--dt",type=float,default=1.0)
    pa.add_argument("--out",type=str,default="roughness_lnsv")
    pa.add_argument("--fast",action="store_true",help="H=[0.05,0.10,0.40], lam=[1.0], n=8, T=800")
    args=pa.parse_args()
    if args.fast:
        args.H=[0.05,0.10,0.40]; args.lam=[1.0]; args.lam_fixed=1.0
        args.n_dgp=8; args.n_cal=3; args.n_sim=8; args.T=800; args.burn=200
    run(H_grid=args.H,lam_grid=args.lam,lam_fixed=args.lam_fixed,
        n_dgp=args.n_dgp,n_sim=args.n_sim,T=args.T,burn=args.burn,
        dt=args.dt,n_cal=args.n_cal,T_cal=min(args.T,600),burn_cal=min(args.burn,150),
        save_prefix=args.out)
