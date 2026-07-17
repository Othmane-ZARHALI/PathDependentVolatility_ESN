"""
zumbach_pdvgl.py
================

DGP: the Guyon & Lekeufack path-dependent volatility model (PDV-GL),
in its 4-FACTOR MARKOVIAN form (power-law kernels approximated by a
mixture of two exponentials each), used HERE AS THE GROUND-TRUTH DGP --
not as a candidate model to be calibrated, as in earlier drafts of this
protocol.

    dS_t/S_t = sigma_t dW_t
    sigma_t  = beta0 + beta1*R1_t + beta2*sqrt(R2_t) + beta12*R1_t*sqrt(R2_t)
    R1_t = theta1*X_{1,0,t} + (1-theta1)*X_{1,1,t}   (leverage factor,
           X_{1,j} = EWMA of past signed returns, rate lam1j)
    R2_t = theta2*X_{2,0,t} + (1-theta2)*X_{2,1,t}   (magnitude factor,
           X_{2,j} = EWMA of past squared returns, rate lam2j)

"True" (ground-truth) parameters: literature-calibrated "Table 8,
Parameter Set 1 (Realistic paths)" --
    beta0=0.04, beta1=-0.13, beta2=0.65, beta12=0 (unused for Set 1)
    lam10=55, lam11=10, theta1=0.25
    lam20=20, lam21=3,  theta2=0.5
Run literally these give ~4% annualised vol (return-scale convention
uncertain from the table alone); a documented post-hoc linear rescale
to the common 20% target is applied (see dgp_pdvgl docstring). This is
a single, fully-specified calibration -- no swept grid parameter.

ONLY THE ESN AND QRH ARE COMPARED against this DGP (PDV-GL-as-a-model is
deliberately NOT run: PDV-GL is now the ground truth itself, so
comparing it against a freshly-calibrated copy of itself would be
close to tautological).

MODEL: the optimal architecture found by the ESN hyperparameter-search
protocol (V5): Nr=64, Nz=8, 9 shared hyperparameters from the smooth-
score grid-matrix sweep. H, lam_lo, lam_hi, az_lo, az_hi are CALIBRATED
per alpha_r cell (not fixed) -- see quick_calibrate.

CALIBRATION: full 11-term data-adaptive score S, with the SMOOTH
(Gaussian-kernel) proximity of the ESN hyperparameter-search protocol V5
(score_fn), replacing the original piecewise-linear tent kernel.
S >= 0, MINIMISED directly, S=0 is a perfect match.

ZUMBACH FIGURES: unchanged in structure (Z1-Z4 + main cross-moment
figure), now with only 2 sources (DGP + ESN) instead of 4; "H"/"eta"
axis labels re-purposed to mean alpha_r (true_alpha_r is fixed per
column, eta_grid is a single dummy value).

OUTPUT FILES (prefix zumbach_pdvgl):
  zumbach_pdvgl_Z1_{w}.png
  zumbach_pdvgl_Z2.png
  zumbach_pdvgl_Z3_{w}.png
  zumbach_pdvgl_Z4.png
  zumbach_pdvgl_zumbach_{w}.png
  zumbach_pdvgl_protocol.tex / .pdf
"""

import math, warnings, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats, optimize
from scipy.special import gamma as Gamma
from scipy.linalg import cholesky as sp_chol
from scipy.ndimage import gaussian_filter1d

warnings.filterwarnings("ignore")

# ============================================================
# 0.  Constants and grid
# ============================================================

TRADING_DAYS = 252.0
TV_ANN  = 0.20
TV_DAY  = TV_ANN / math.sqrt(TRADING_DAYS)
SIG_MIN = 0.01 / math.sqrt(TRADING_DAYS)

TRUE_ALPHA_R_GRID = [1]   # dummy single-value "grid": Parameter Set 1 is
                          # one fixed literature calibration, no swept parameter
DGP_WEIGHT_FIXED  = None            # placeholder, unused (kept for CLI compat)
# kept for backward-compatible variable names used throughout run()/CLI:
H_GRID, ETA_GRID = TRUE_ALPHA_R_GRID, [None]
RHO_DGP  = None           # not used; PDV-GL's own leverage comes from its
                          # beta1*R1 (signed past returns) term, not a
                          # separate rho parameter

# "True" (ground-truth) PDV-GL parameters used to generate the DGP itself.
# alpha_r (the R1/leverage-kernel exponent) is swept over
# TRUE_ALPHA_R_GRID; the others are held fixed. beta1 is pushed to a
# STRONG, clearly-negative value so the DGP shows a strong, sustained
# Zumbach/time-reversal-asymmetry effect (matching the requested
# "strong Zumbach effect" reference), rather than the barely-there effect
# of a mild leverage coefficient.
TRUE_BETA0   = TV_DAY * 0.75   # floor/baseline daily vol level (tuned so
                               # the realised annualised vol lands near
                               # the 20% target once the feedback terms
                               # are added on top)
TRUE_BETA1   = -1.20    # leverage (R1) loading -- NEGATIVE => genuine
                        # leverage/Zumbach asymmetry (past negative
                        # return raises future vol); tuned for a STRONG,
                        # sustained cross-moment Zumbach effect
                        # (|Z(tau)/sigma^3| ~ 0.03-0.09 for tau=1..30,
                        # comparable in magnitude to the reference
                        # "PDV-GL (Table 8, Set 1)" bars)
TRUE_BETA2   = 0.28     # magnitude (R2, sqrt) loading
TRUE_ALPHA_V = 0.25     # long-memory (R2) kernel exponent, fixed

COLORS = {
    "DGP (PDV-GL)":   "#888780",
    "ESN A2-981003":  "#1a4d80",
    "QRH":            "#D85A30",
}
CORAL = "#c0392b"; GRAY = "#888780"
MODEL_NAMES = ["ESN A2-981003", "QRH"]   # ESN and QRH compared against the DGP

_ARCH = dict(
    matrix_seed=202695547565, n_r=96, n_z=12, H_target=0.08,
    z_strength=0.19288535085059877, even_strength=2.2291053209265703,
    linear_strength=0.1930215632408165, gamma_norm=0.9651266728165394,
    local_z_strength=0.04759055003377607, zz_scale=0.039730814766135034,
    sign_prob_neg=0.22189311131648579, rough_orientation=-1.0,
)
# _ARCH = the ESN model's current optimal architecture: Nr=96, Nz=12,
# with 7 shared hyperparameters. H_target is a start-up default only;
# H, lam_lo, lam_hi, az_lo, az_hi, rough_scale, zr_lo, zr_hi, m1, m2 are
# all CALIBRATED per (alpha,weight) DGP cell (quick_calibrate), not
# fixed. The read-out is
#   eta_t = b0 + rho_r*rough_scale*(q^T r_t) + sum_j w_{j,z} z_{j,t}
#           + r_t^T Q r_t,
#   w_{j,z} = z_readout_j/sqrt(Nz), z_readout_j = geomspace(zr_lo,zr_hi,Nz)[j],
#   Q = (1/Nr)*(m1*I + m2*q q^T)  (symmetric, no PD constraint).
ROUGH_SCALE_BOUNDS = (0.05, 2.0)
ZR_LO_BOUNDS = (0.01, 0.5)
ZR_HI_BOUNDS = (0.01, 0.5)
M1_BOUNDS = (-2.0, 2.0)
M2_BOUNDS = (-2.0, 2.0)

# ============================================================
# 1.  Shared statistics (IDENTICAL to roughness_evaluation.py)
# ============================================================

def _sp(x):
    if x>35: return x
    if x<-35: return math.exp(x)
    return math.log1p(math.exp(x))

def _inv_sp(y):
    y=max(float(y),1e-15)
    return y if y>35 else math.log(math.expm1(y))

def _corr(a, b):
    a,b=np.asarray(a,float),np.asarray(b,float)
    if len(a)<3 or a.std()<1e-14 or b.std()<1e-14: return np.nan
    return float(np.corrcoef(a,b)[0,1])

def compute_statistics(daily_x, daily_var):
    """11 stylized-fact statistics — identical to roughness_evaluation.py."""
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
    zvals=[]; t1vals=[]; t2vals=[]
    for L in [5,10,20]:
        n=len(x)-2*L
        if n>30:
            csx=np.concatenate([[0.],np.cumsum(x)])
            csv=np.concatenate([[0.],np.cumsum(v)])
            idx=L+np.arange(n)
            pR=csx[idx]-csx[idx-L]; fR=csx[idx+L]-csx[idx]
            pV=(csv[idx]-csv[idx-L])/L; fV=(csv[idx+L]-csv[idx])/L
            t1=_corr(pR**2,fV); t2=_corr(pV,fR**2)
            t1vals.append(t1); t2vals.append(t2); zvals.append(t1-t2)
    c=x-x.mean(); var=float(np.var(c))
    kurt=float(np.mean(c**4)/(var**2+1e-30)-3.)
    ann=math.sqrt(TRADING_DAYS)
    return dict(
        H_hat=float(H), mean_vol_ann=float(sigma_daily.mean())*ann,
        q995_vol_ann=float(np.quantile(sigma_daily,0.995))*ann,
        max_vol_ann=float(sigma_daily.max())*ann,
        mean_vol_acf=float(np.nanmean(vol_acf)),
        taylor_gap=float(np.nanmean(abs_acf)-np.nanmean(sq_acf)),
        taylor_frac=float(np.nanmean(abs_acf>sq_acf)),
        zumbach=float(np.nanmean(zvals)),
        zumbach_term1=float(np.nanmean(t1vals)),   # Corr(pastR^2, futV)
        zumbach_term2=float(np.nanmean(t2vals)),   # Corr(pastV, futR^2)
        leverage=float(np.nanmean(lev)), kurtosis=kurt,
        max_ret_acf=float(np.nanmax(np.abs(ret_acf))),
        logret=x, sigma_daily=sigma_daily,
        vol_acf_vec=vol_acf, lev_vec=lev,
    )

def make_score_ref(dgp_stats):
    """Data-adaptive score reference — identical to roughness_evaluation.py."""
    def _m(k): return float(np.nanmean([s[k] for s in dgp_stats]))
    def _sd(k): return float(np.nanstd( [s[k] for s in dgp_stats]))
    floors=dict(H_hat=0.01,mean_vol_ann=0.01,q995_vol_ann=0.05,
                max_vol_ann=0.10,mean_vol_acf=0.005,taylor_gap=0.002,
                taylor_frac=0.05,zumbach=0.005,
                zumbach_term1=0.02,zumbach_term2=0.02,
                leverage=0.005,
                kurtosis=0.20,max_ret_acf=0.005)
    ref={}
    for k,fl in floors.items():
        c=_m(k); s=_sd(k)
        ref[k+"_c"]=c; ref[k+"_s"]=max(s,abs(c)*0.30,fl)
    ref["stress_mx"]=max(ref["max_vol_ann_c"]*2.0,1.5)
    return ref

# ZUMBACH-FIT FIX (part 2): the score originally only targeted the
# SCALAR DIFFERENCE Z=Term1-Term2, which a model can match by accident
# even when BOTH individual channels are wrong -- e.g. ESN and QRH were
# both found to have Term1 (Corr(pastR^2,futV)) AND Term2
# (Corr(pastV,futR^2)) elevated well above the DGP's own values
# (whose Term2 sits near 0, i.e. the DGP has essentially NO
# past-vol-predicts-future-magnitude channel at all), yet their
# DIFFERENCE still landed close to the DGP's Z because both channels
# were inflated together. To close this, the single zumbach-difference
# term is replaced by TWO separate terms, each individually pulling the
# model's Term1 and Term2 toward the DGP's own Term1_c, Term2_c -- a
# model can no longer "cheat" by getting the difference right for the
# wrong reason. The combined weight (3.5, split 1.75/1.75) is unchanged
# from the single zumbach-difference weight.
W_ZTERM1 = 1.75   # Corr(past R^2, fut V)  -- was folded into W_ZUMBACH=3.5
W_ZTERM2 = 1.75   # Corr(past V, fut R^2)  -- was folded into W_ZUMBACH=3.5
SCORE_MAX = 2.2+1.4+0.6+0.2+0.8+1.0+0.8+(W_ZTERM1+W_ZTERM2)+1.1+0.7+0.5 + 5.0
# = 17.8, worst possible S (unchanged: same total Zumbach-related weight
# as the single-term version, just split across the two channels)

def score_fn(st, ref):
    """
    Full 11-term score S, SMOOTH (Gaussian-kernel) version -- EXACT
    convention of the ESN hyperparameter-search protocol (V5): S >= 0,
    MINIMISED, S = 0 <=> a perfect match on every stylised fact and no
    stress event. Replaces the piecewise-linear tent kernel (which hits
    exactly 0, zero gradient, beyond one tolerance width) with
        f_smooth(x,c,s) = exp(-0.5*((x-c)/s)^2)   in (0,1]
        g_smooth(x,c,s) = exp(-0.5*(max(0,x-c)/s)^2)
    which never flattens, giving Nelder-Mead a usable gradient at any
    distance from target.

    MODIFIED (Zumbach-fit fix, part 2): the single zumbach-DIFFERENCE
    term (weight 3.5) is replaced by TWO terms matching Term1
    (Corr(pastR^2,futV), weight 1.75) and Term2 (Corr(pastV,futR^2),
    weight 1.75) INDIVIDUALLY -- see the rationale directly above.
    """
    def f(x,c,s): return math.exp(-0.5*((x-c)/max(s,1e-8))**2)
    def g(x,c,s): return math.exp(-0.5*(max(0.,x-c)/max(s,1e-8))**2)
    H=st["H_hat"]; vol=st["mean_vol_ann"]; q995=st["q995_vol_ann"]; mx=st["max_vol_ann"]
    V=st["mean_vol_acf"]; GT=st["taylor_gap"]; FT=st["taylor_frac"]
    Z1=st["zumbach_term1"]; Z2=st["zumbach_term2"]
    L=st["leverage"]; K=st["kurtosis"]; A=st["max_ret_acf"]
    stress=int(mx>ref["stress_mx"] or vol<0.05 or vol>1.50)
    d =2.2*(1.-f(H,  ref["H_hat_c"],       ref["H_hat_s"]))
    d+=1.4*(1.-f(vol,ref["mean_vol_ann_c"], ref["mean_vol_ann_s"]))
    d+=0.6*(1.-g(q995,ref["q995_vol_ann_c"],ref["q995_vol_ann_s"]))
    d+=0.2*(1.-g(mx, ref["max_vol_ann_c"],  ref["max_vol_ann_s"]))
    d+=0.8*(1.-f(V,  ref["mean_vol_acf_c"], ref["mean_vol_acf_s"]))
    d+=1.0*(1.-f(GT, ref["taylor_gap_c"],   ref["taylor_gap_s"]))
    d+=0.8*(1.-f(FT, ref["taylor_frac_c"],  ref["taylor_frac_s"]))
    d+=W_ZTERM1*(1.-f(Z1, ref["zumbach_term1_c"], ref["zumbach_term1_s"]))
    d+=W_ZTERM2*(1.-f(Z2, ref["zumbach_term2_c"], ref["zumbach_term2_s"]))
    d+=1.1*(1.-f(L,  ref["leverage_c"],     ref["leverage_s"]))
    d+=0.7*(1.-f(K,  ref["kurtosis_c"],     ref["kurtosis_s"]))
    d+=0.5*(1.-f(A,  ref["max_ret_acf_c"],  ref["max_ret_acf_s"]))
    d+=5.0*stress
    return float(d)

def hill_curve(x, n_pts=40):
    n=len(x); sx=np.sort(np.abs(x))[::-1]
    ks=np.unique(np.round(np.exp(
        np.linspace(np.log(5),np.log(max(n//5,6)),n_pts))).astype(int))
    out=[]
    for k in ks:
        k=min(k,n-1); lr=np.log(sx[:k])-np.log(sx[k]); lr=lr[lr>0]
        out.append(1./np.mean(lr) if len(lr)>0 else np.nan)
    return ks, np.array(out)

# ============================================================
# 2.  DGP — PDV-GL (ground truth)
# ============================================================

_RL_CACHE = {}

def _get_rl_chol(T, H, dt=1.0):
    """
    Cholesky factor of the Riemann-Liouville covariance matrix.

    The RL process W^H_t = sqrt(2H) * int_0^t (t-s)^{H-1/2} dW_s has
    covariance:
        Cov(W^H_s, W^H_t) = 2H * int_0^{min(s,t)} (s-u)^{H-0.5}(t-u)^{H-0.5} du

    For simulation efficiency we use the fBM covariance as a proxy:
        Sigma_ij = 0.5*(t_i^{2H} + t_j^{2H} - |t_i-t_j|^{2H})
    which has the same Holder exponent H and identical short-lag variance
    structure (both have Var(increment over [s,s+tau]) ~ tau^{2H}).
    The Cholesky factor is cached per (T, H, dt).
    """
    key=(T,H,dt)
    if key not in _RL_CACHE:
        t=np.arange(1,T+1,dtype=float)*dt
        C=0.5*(t[:,None]**(2*H)+t[None,:]**(2*H)
               -np.abs(t[:,None]-t[None,:])**(2*H))
        C+=np.eye(T)*1e-10
        try:    L=sp_chol(C,lower=True)
        except: L=np.linalg.cholesky(C+np.eye(T)*1e-8)
        _RL_CACHE[key]=L
    return _RL_CACHE[key]

# ============================================================
# 3.  Models — UNCHANGED from roughness_evaluation.py
# ============================================================

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
    """Build reservoir matrices. H, lam_lo, lam_hi, az_lo, az_hi,
    rough_scale, zr_lo, zr_hi, m1, m2 are the per-(alpha,weight)-cell
    CALIBRATED inner-loop parameters; if H is omitted, falls back to
    arch['H_target'] (start-up admissibility check only).
    The read-out generalises z_readout to a per-mode profile
    z_readout_j = geomspace(zr_lo, zr_hi, Nz)[j], and adds a quadratic
    term r_t^T Q r_t, Q=(1/Nr)(m1*I+m2*qq^T) -- a genuine SYMMETRIC
    matrix (no positive-definiteness constraint) -- to the
    pre-activation (see _sim_esn)."""
    if H is None: H = arch["H_target"]
    rng=np.random.default_rng(int(arch["matrix_seed"]))
    n_r,n_z=arch["n_r"],arch["n_z"]
    lam,q=_kernel_nodes(n_r,H,lam_lo,lam_hi)
    az_lo=max(float(az_lo),1e-8); az_hi=max(float(az_hi),az_lo*1.5)
    az=np.geomspace(az_lo,az_hi,n_z)
    zr_lo_c=max(float(zr_lo),1e-4); zr_hi_c=max(float(zr_hi),zr_lo_c*1.05)
    z_readout_profile=np.geomspace(zr_lo_c,zr_hi_c,n_z)
    zz=rng.uniform(-arch["zz_scale"],arch["zz_scale"],n_z)
    sz=-rng.choice([-1.,1.],n_z,p=[arch["sign_prob_neg"],1-arch["sign_prob_neg"]])
    fi=np.array([min(n_r-1,n_r//2+int((n_r//2-1)*j/max(n_z-1,1))) for j in range(n_z)],dtype=int)
    si=np.array([int((n_r//2-1)*(n_z-1-j)/max(n_z-1,1)) for j in range(n_z)],dtype=int)
    b0=_inv_sp(math.sqrt(max(TV_DAY**2-SIG_MIN**2,1e-15)))
    k0=arch["rough_orientation"]*rough_scale*float(q@np.sqrt(2.*lam))
    return dict(lam=lam,q=q,az=az,zz=zz,sz=sz,fi=fi,si=si,b0=b0,k0=k0,
                rough_scale=float(rough_scale),
                z_readout_profile=z_readout_profile, m1=float(m1), m2=float(m2))

_ESN_PARAMS=None

def _sim_esn(seed,T,dt=1.0,arch=_ARCH,params=None):
    """eta_t = b0 + rho_r*rough_scale*(q^T r_t) + sum_j w_{j,z}*z_{j,t}
    + r_t^T Q r_t, Q=(1/Nr)*(m1*I+m2*qq^T), w_{j,z} =
    z_readout_j/sqrt(Nz). rough_scale, the z_readout_profile, and
    (m1,m2) are read from `params` (calibrated per DGP cell)."""
    global _ESN_PARAMS
    if params is None:
        if _ESN_PARAMS is None: _ESN_PARAMS=_build_esn(arch)
        params=_ESN_PARAMS
    P=params; rng=np.random.default_rng(int(seed))
    n_r,n_z=arch["n_r"],arch["n_z"]
    spd=int(round(1./dt)); n_st=T*spd; sdt=math.sqrt(dt)
    al=np.exp(-P["lam"]*dt); cl=np.sqrt(np.maximum(1.-al**2,1e-14))
    azd=np.exp(-P["az"]*dt); om=1.-azd
    r=np.zeros(n_r); z=np.zeros(n_z); dx=np.zeros(T); dv=np.zeros(T)
    wz_vec=P["z_readout_profile"]/math.sqrt(n_z)
    rc=arch["rough_orientation"]*P["rough_scale"]
    m1c=P["m1"]/n_r; m2c=P["m2"]/n_r
    qvec=P["q"]
    for step in range(n_st):
        eps=rng.normal()
        qr=float(qvec@r)
        quad_term=m1c*float(r@r)+m2c*qr*qr
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

def _sim_gl(seed,T,dt,beta0,beta1,beta2,alpha_r,alpha_v):
    """
    CORRECTED Guyon & Lekeufack (2023) two-factor model:
        sigma_t = beta0 + beta1*R1_t + beta2*sqrt(R2_t)
        R1_t = power-law-weighted sum of PAST RAW (signed) returns
               -- kernel exponent alpha_r; this is the LEVERAGE factor:
               beta1 < 0 means a past NEGATIVE return raises future vol.
        R2_t = power-law-weighted sum of PAST SQUARED returns
               -- kernel exponent alpha_v; this is the MAGNITUDE
               (vol-clustering) factor.
    An earlier draft of this DGP used two power-law-weighted-SQUARED-
    return factors (no signed R1 term at all), which is symmetric in the
    sign of returns and therefore CANNOT produce a non-zero cross-moment
    Zumbach statistic Z^cm(tau)=E[r^2_{t+tau}r_t]-E[r_{t+tau}r_t^2] for
    ANY parameter choice (both factors are even functions of r). This
    corrected version restores the genuine leverage/asymmetry term
    (R1, linear in signed past returns) that the actual Guyon-Lekeufack
    model has, which is what produces a real, non-trivial Zumbach effect.
    """
    spd=int(round(1./dt)); n_st=T*spd; sdt=math.sqrt(dt)
    rng=np.random.default_rng(int(seed))
    dx=np.zeros(T); dv=np.zeros(T); rbuf=[]
    sig=beta0
    for step in range(n_st):
        eps=rng.normal()
        r_t=sig*sdt*eps; day=step//spd
        dx[day]+=r_t; dv[day]+=sig*sig*dt
        rbuf.append(r_t); t=len(rbuf)
        if t<2:
            sig=beta0
        else:
            ra=np.array(rbuf); lgs=np.arange(1,t+1,dtype=float)
            r_rev=ra[::-1]/sdt                # raw (unit-time) past returns
            r2_rev=r_rev**2
            K1=lgs**(alpha_r-0.5); R1=float((K1*r_rev).sum()/K1.sum())
            K2=lgs**(alpha_v-0.5); R2=float((K2*r2_rev).sum()/K2.sum())
            sig=max(beta0+beta1*R1+beta2*math.sqrt(max(R2,0.)), 1e-5)
    return dx,dv

# ── 4-factor Markovian PDV-GL (Guyon & Lekeufack, "Table 8, Parameter
#    Set 1 (Realistic paths)" calibration) ──────────────────────────────
# Each of the power-law R1 (leverage) and R2 (magnitude) factors of
# _sim_gl above is here replaced by a MARKOVIAN approximation: a mixture
# of two exponential (OU-type) kernels, i.e. a "4-factor" model (2
# factors for R1, 2 for R2). This is the standard trick used to make
# genuinely path-dependent power-law-kernel volatility models fast to
# simulate as a finite-dimensional Markov process.
DT_YEARS = 1.0 / TRADING_DAYS

def _sim_gl_4factor(seed, T, beta0, beta1, beta2, beta12,
                     lam10, lam11, theta1, lam20, lam21, theta2, burn=0):
    """
    sigma_t = beta0 + beta1*R1_t + beta2*sqrt(R2_t) + beta12*R1_t*sqrt(R2_t)
    R1_t = theta1*X_{1,0,t} + (1-theta1)*X_{1,1,t}      (leverage factor)
    R2_t = theta2*X_{2,0,t} + (1-theta2)*X_{2,1,t}      (magnitude factor)
    where X_{1,j,t} is an exponential (EWMA) factor of past SIGNED
    returns with rate lam1j, and X_{2,j,t} is an exponential factor of
    past SQUARED returns with rate lam2j (both rates in ANNUALISED
    units, matched to the literature calibration's own convention;
    decayed at each daily step via exp(-lam*DT_YEARS)).
    Returns (r, v) raw arrays of length T+burn; caller discards burn.
    """
    rng = np.random.default_rng(int(seed))
    N = T + burn
    x10=x11=x20=x21=0.0
    e10=math.exp(-lam10*DT_YEARS); e11=math.exp(-lam11*DT_YEARS)
    e20=math.exp(-lam20*DT_YEARS); e21=math.exp(-lam21*DT_YEARS)
    r=np.zeros(N); v=np.zeros(N)
    for t in range(N):
        R1 = theta1*x10 + (1-theta1)*x11
        R2 = theta2*x20 + (1-theta2)*x21
        sqR2 = math.sqrt(max(R2, 0.))
        sig = max(beta0 + beta1*R1 + beta2*sqR2 + beta12*R1*sqR2, 1e-6)
        eps = rng.standard_normal()
        rt = sig * math.sqrt(DT_YEARS) * eps
        r[t] = rt; v[t] = sig*sig*DT_YEARS
        x10 = e10*x10 + (1-e10)*rt
        x11 = e11*x11 + (1-e11)*rt
        x20 = e20*x20 + (1-e20)*rt*rt
        x21 = e21*x21 + (1-e21)*rt*rt
    return r, v

# "True" (ground-truth) parameters: literature-calibrated Table 8,
# "Parameter Set 1 (Realistic paths)" of the 4-factor Markovian PDV
# model. Run through _sim_gl_4factor literally, the realised annualised
# vol comes out at only ~4% (the table's implicit return-scale
# convention could not be verified from the parameters alone), so a
# single documented linear RESCALING by TV_DAY/realised_std is applied
# post-hoc to bring the level to the same 20% annualised target used by
# every DGP in this line of protocols -- this rescaling is a pure
# amplitude change and leaves every normalised statistic (Hurst,
# kurtosis, the Z(tau)/sigma^3 cross-moment ratio, leverage correlation)
# EXACTLY unchanged, since numerator and denominator scale identically.
TRUE_4F_BETA0  = 0.04
TRUE_4F_BETA1  = -0.13
TRUE_4F_BETA2  = 0.65
TRUE_4F_BETA12 = 0.0     # unused for Set 1
TRUE_4F_LAM10  = 55.0
TRUE_4F_LAM11  = 10.0
TRUE_4F_THETA1 = 0.25
TRUE_4F_LAM20  = 20.0
TRUE_4F_LAM21  = 3.0
TRUE_4F_THETA2 = 0.5

def dgp_pdvgl(alpha_r, n_paths, T, burn, dt=1.0, seed_base=0, **kwargs):
    """
    Guyon & Lekeufack path-dependent volatility model, 4-FACTOR MARKOVIAN
    form with the literature "Table 8, Parameter Set 1" calibration, used
    HERE AS THE DGP ITSELF (ground truth), not as a model to calibrate.
    `alpha_r` is accepted purely for interface compatibility with the
    rest of this pipeline (single-column loop variable) -- Parameter
    Set 1 has no free "alpha_r" of its own; every path uses the SAME
    fixed literature parameters (TRUE_4F_*).
    """
    rng_base = np.random.default_rng(seed_base)
    st_all = []
    for p in range(n_paths):
        seed = int(rng_base.integers(1 << 31))
        r_raw, v_raw = _sim_gl_4factor(
            seed, T, TRUE_4F_BETA0, TRUE_4F_BETA1, TRUE_4F_BETA2, TRUE_4F_BETA12,
            TRUE_4F_LAM10, TRUE_4F_LAM11, TRUE_4F_THETA1,
            TRUE_4F_LAM20, TRUE_4F_LAM21, TRUE_4F_THETA2, burn=burn)
        r_used = r_raw[burn:]; v_used = v_raw[burn:]
        realised_std = r_used.std()
        scale = (TV_DAY / realised_std) if realised_std > 1e-14 else 1.0
        X = r_used * scale
        V = v_used * scale**2
        st = compute_statistics(X, V)
        st_all.append(st)
    return st_all

def sim_model_paths(model_name,cal,n_paths,T,burn,dt=1.0):
    """Simulate one model and compute all 11 score statistics."""
    st_all=[]
    esn_params=None
    if model_name=="ESN A2-981003":
        esn_params=_build_esn(_ARCH, H=cal["H"], lam_lo=cal["lam_lo"], lam_hi=cal["lam_hi"],
                               az_lo=cal["az_lo"], az_hi=cal["az_hi"],
                               rough_scale=cal.get("rough_scale",0.40),
                               zr_lo=cal.get("zr_lo",0.03), zr_hi=cal.get("zr_hi",0.07),
                               m1=cal.get("m1",0.0), m2=cal.get("m2",0.0))
        esn_params["b0"]=esn_params["b0"]+cal.get("b0_delta",0.)
    for p in range(n_paths):
        if model_name=="ESN A2-981003":
            x,v=_sim_esn(9000+p,T,dt,arch=_ARCH,params=esn_params); x=x*cal.get("scale",1.0)
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
    Calibrate the ESN and QRH to the PDV-GL DGP using the FULL 11-term
    smooth (Gaussian-kernel) data-adaptive score S. ESN calibrates its
    full 12-dim inner-loop vector (H, lam_lo, lam_hi, az_lo, az_hi,
    b0_delta, scale, rough_scale, zr_lo, zr_hi, m1, m2) via 11-start
    Nelder-Mead; only the 7 shared architecture hyperparameters
    (_ARCH) are held fixed.
    PDV-GL-as-a-calibrated-model is DELIBERATELY NOT run here: PDV-GL is
    now the ground-truth DGP itself, so comparing it against a freshly
    -fit copy of itself would add little (see run()).
    """
    ref=make_score_ref(dgp_sts); tH=ref["H_hat_c"]
    print(f"    Score ref: H={tH:.4f}  vol={ref['mean_vol_ann_c']*100:.1f}%  "
          f"Z={ref['zumbach_c']:.4f}  L={ref['leverage_c']:.4f}  K={ref['kurtosis_c']:.3f}")

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
            params=_build_esn(_ARCH, H=H, lam_lo=lam_lo, lam_hi=lam_hi, az_lo=az_lo, az_hi=az_hi,
                               rough_scale=rough_scale, zr_lo=zr_lo, zr_hi=zr_hi,
                               m1=m1, m2=m2)
            if params["k0"] >= 0: return SCORE_MAX
            params["b0"]=params["b0"]+b0d
            scores=[]
            for q in range(n_cal):
                x,v=_sim_esn(3000+q,T_cal,dt,arch=_ARCH,params=params)
                xb=(x*sc)[burn_cal:]; vb=v[burn_cal:]
                scores.append(score_fn(compute_statistics(xb,vb),ref))
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
        (0.10, 1e-5,   5.0,   0.3,   0.5,  0.40, 0.03, 0.07, 0.0, 0.0),  # Zumbach term1/term2-balanced
        (0.15, 1e-4,   1.0,   1/300, 1/8,  0.90, 0.10, 0.30, 0.0, 0.0),  # strong amplitude
        (0.20, 1e-4,   1.0,   1/300, 1/8,  0.15, 0.02, 0.03, 0.0, 0.0),  # weak amplitude
        (0.20, 1e-4,   1.0,   1/300, 1/8,  0.40, 0.02, 0.20, 0.8, 0.0),  # isotropic +
        (0.20, 1e-4,   1.0,   1/300, 1/8,  0.40, 0.20, 0.02, -0.8, 0.0), # isotropic -
        (0.20, 1e-4,   1.0,   1/300, 1/8,  0.40, 0.03, 0.07, 0.0, 0.8),  # q-aligned +
    ]
    best_fun, best_x = None, None
    n_starts = len(ESN_INNER_STARTS)
    maxiter_per_start = max(30, 330 // n_starts)
    for (H0, ll0, lh0, al0, ah0, rs0, zrl0, zrh0, m10, m20) in ESN_INNER_STARTS:
        x0=[H0, math.log(ll0), math.log(lh0), math.log(al0), math.log(ah0), 0.0, 0.0,
            math.log(rs0), math.log(zrl0), math.log(zrh0), m10, m20]
        res=optimize.minimize(esn_obj, x0, method="Nelder-Mead",
                              options={"maxiter":maxiter_per_start,"xatol":0.02,"fatol":0.02})
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
    esn_cal=dict(H=H, lam_lo=lam_lo, lam_hi=lam_hi, az_lo=az_lo, az_hi=az_hi,
                 b0_delta=b0d, scale=sc, rough_scale=rough_scale,
                 zr_lo=zr_lo, zr_hi=zr_hi, m1=m1, m2=m2)

    def qrh_obj(p):
        H,nu,lk,cf=p
        if not(0.01<H<0.49 and 0.01<nu<3. and 0.1<lk<20. and 0.01<cf<0.99): return SCORE_MAX
        scores=[]
        for q in range(n_cal):
            try:
                x,v=_sim_qrh(500+q,T_cal,dt,H,nu,lk,cf)
                scores.append(score_fn(compute_statistics(x[burn_cal:],v[burn_cal:]),ref))
            except: scores.append(SCORE_MAX)
        return float(np.nanmean(scores))
    # MODIFIED (Zumbach-fit fix): multi-start (was a single start at
    # [tH,0.30,2.0,0.50], which is biased toward the DGP's H_hat target
    # and gets stuck there). Sensitivity testing showed QRH's Zumbach
    # effect is strongly anti-correlated with H: low H (~0.03-0.10)
    # gives a robust, DGP-matching positive Zumbach (0.14-0.24), so a
    # dedicated low-H start is added alongside the original H_hat-biased
    # start, keeping the best-of-both result.
    QRH_STARTS = [
        [max(0.02,min(tH,0.48)), 0.30, 2.0, 0.50],   # original: H_hat-biased
        [0.06, 0.40, 1.5, 0.40],                      # low-H: Zumbach-favorable
        [0.06, 0.15, 0.5, 0.80],                      # low-H, low-nu_vol: balances
                                                        # Term1 AND Term2 individually
                                                        # (see zumbach-results section)
    ]
    best_fun, best_x = None, None
    for x0 in QRH_STARTS:
        res=optimize.minimize(qrh_obj,x0,method="Nelder-Mead",
                              options={"maxiter":300,"xatol":0.01,"fatol":0.02})
        if best_fun is None or res.fun < best_fun:
            best_fun, best_x = res.fun, res.x
    H_o,nv_o,lk_o,cf_o=best_x
    H_o=float(np.clip(H_o,0.01,0.48)); nv_o=float(np.clip(nv_o,0.01,3.))
    lk_o=float(np.clip(lk_o,0.1,20.)); cf_o=float(np.clip(cf_o,0.01,0.99))
    print(f"    QRH: H={H_o:.4f}  nu={nv_o:.4f}  lam={lk_o:.4f}  c={cf_o:.4f}  score={best_fun:.3f}")
    qrh_cal=dict(H=H_o,nu_vol=nv_o,lam=lk_o,c_frac=cf_o)

    return esn_cal, qrh_cal, ref

# ============================================================
# 5.  Zumbach utility: per-lag Z(L) curve
# ============================================================

def zumbach_curve(st_list, lags=None):
    """
    Z(L) = Corr(past_R^2, future_V) - Corr(past_V, future_R^2)
    computed for each lag L and each path, then averaged cross-path.

    Returns: Zm (mean), Zs (std), shape (len(lags),)
    """
    if lags is None: lags=list(range(1,31))
    all_curves=[]
    for st in st_list:
        x=st["logret"]; v=st["sigma_daily"]**2
        curve=[]
        for L in lags:
            n=len(x)-2*L
            if n>30:
                csx=np.concatenate([[0.],np.cumsum(x)])
                csv=np.concatenate([[0.],np.cumsum(v)])
                idx=L+np.arange(n)
                pR=csx[idx]-csx[idx-L]; fR=csx[idx+L]-csx[idx]
                pV=(csv[idx]-csv[idx-L])/L; fV=(csv[idx+L]-csv[idx])/L
                curve.append(_corr(pR**2,fV)-_corr(pV,fR**2))
            else:
                curve.append(np.nan)
        all_curves.append(curve)
    arr=np.array(all_curves,dtype=float)
    return np.nanmean(arr,axis=0), np.nanstd(arr,axis=0)


def zumbach_windowed_decomp(st_list, max_lag=30):
    """
    Per-path, per-lag WINDOWED Zumbach decomposition (the "good" metric):
        Z(L)  = Corr(past_R^2, future_V) - Corr(past_V, future_R^2)
        C1(L) = Corr(past_R^2, future_V)   (past-squared-return -> future-vol
                 channel; large POSITIVE values = genuine Zumbach direction)
        C2(L) = Corr(past_V, future_R^2)   (past-vol -> future-squared-return
                 channel; the "reverse-time" direction)
    Z(L) > 0 <=> the past-magnitude -> future-vol channel dominates the
    reverse channel = Zumbach / time-reversal-asymmetry effect present.
    (This replaces the cross-moment Z^cm(tau)=E[r^2_{t+tau}r_t]-E[r_{t+tau}r_t^2],
    which requires a signed leverage effect and is a DIFFERENT, less
    universal diagnostic than the windowed correlation-based Z(L) used for
    calibration -- see score_fn / compute_statistics.)

    Returns: Zm, Zs, C1m, C1s, C2m, C2s, each shape (max_lag,)
    """
    lags = np.arange(1, max_lag + 1)
    Z_all=[]; C1_all=[]; C2_all=[]
    for st in st_list:
        x = st["logret"]; v = st["sigma_daily"]**2
        Zc=[]; C1c=[]; C2c=[]
        for L in lags:
            n = len(x) - 2*L
            if n > 30:
                csx = np.concatenate([[0.], np.cumsum(x)])
                csv = np.concatenate([[0.], np.cumsum(v)])
                idx = L + np.arange(n)
                pR = csx[idx]-csx[idx-L]; fR = csx[idx+L]-csx[idx]
                pV = (csv[idx]-csv[idx-L])/L; fV = (csv[idx+L]-csv[idx])/L
                c1 = _corr(pR**2, fV); c2 = _corr(pV, fR**2)
                C1c.append(c1); C2c.append(c2); Zc.append(c1-c2)
            else:
                C1c.append(np.nan); C2c.append(np.nan); Zc.append(np.nan)
        Z_all.append(Zc); C1_all.append(C1c); C2_all.append(C2c)
    Z_a=np.array(Z_all,dtype=float); C1_a=np.array(C1_all,dtype=float); C2_a=np.array(C2_all,dtype=float)
    return (np.nanmean(Z_a,0), np.nanstd(Z_a,0),
            np.nanmean(C1_a,0), np.nanstd(C1_a,0),
            np.nanmean(C2_a,0), np.nanstd(C2_a,0))

# ============================================================
# 6.  Figures
# ============================================================

ALL_NAMES=["DGP (PDV-GL)"]+MODEL_NAMES
SIGMA_SMOOTH=2.0   # Gaussian smoothing sigma for Z(L) curves

def _col(name): return COLORS.get(name,GRAY)

# ── Figure Z1: Z(L) curves, one column per H ─────────────────

def plot_Z1_curves(results, H_grid, eta_fixed, save_path=None):
    """
    Figure Z1 — Zumbach curve Z(L) vs lag L.
    One column per H value (fixed eta).
    One row per source (DGP + 3 models).
    Each panel: per-path mean (smoothed line) + ±1std shaded band + raw dots.
    """
    Hcols=[H for H in H_grid if (H,eta_fixed) in results]
    n_c=len(Hcols); n_r=len(ALL_NAMES)
    lags=list(range(1,31))

    fig,axes=plt.subplots(n_r,n_c,figsize=(5*n_c,4*n_r))
    if n_c==1: axes=axes[:,None]
    if n_r==1: axes=axes[None,:]

    fig.suptitle(
        r"Zumbach Effect $Z(L) = \mathrm{Corr}(\mathrm{past\,}R^2,\,\mathrm{fut\,}V)"
        r"- \mathrm{Corr}(\mathrm{past\,}V,\,\mathrm{fut\,}R^2)$"
        f"\nPDV-GL DGP (ground truth)"
        r"  |  $Z(L)>0$ = Zumbach effect present",
        fontsize=11,fontweight="bold")

    for ci,H in enumerate(Hcols):
        data=results[(H,eta_fixed)]
        axes[0,ci].set_title(f"$\\alpha={H}$",fontsize=11,fontweight="bold")
        for ri,name in enumerate(ALL_NAMES):
            ax=axes[ri,ci]
            Zm,Zs=zumbach_curve(data[name],lags)
            Zm_sm=gaussian_filter1d(np.where(np.isfinite(Zm),Zm,0.),SIGMA_SMOOTH)
            ax.fill_between(lags,Zm-Zs,Zm+Zs,alpha=0.15,color=_col(name))
            ax.scatter(lags,Zm,color=_col(name),s=8,alpha=0.45,zorder=3)
            ax.plot(lags,Zm_sm,color=_col(name),lw=2.,zorder=4,label=name)
            ax.axhline(0,color=GRAY,lw=0.8,ls="--",alpha=0.6)
            ax.set_xlabel("Lag $L$ (days)" if ri==n_r-1 else "",fontsize=9)
            if ci==0: ax.set_ylabel(f"{name}\n$Z(L)$",fontsize=8)
            ax.set_title("" if ri>0 else f"$\\alpha={H}$",fontsize=10 if ri==0 else 9)
            ax.legend(fontsize=7); ax.grid(True,alpha=0.18)

    plt.tight_layout(rect=[0,0,1,0.94])
    if save_path:
        fig.savefig(save_path,dpi=150,bbox_inches="tight")
        print(f"  Saved -> {save_path}")
    return fig


# ── Figure Z2: Zumbach scalar bar chart across H and eta ─────

def plot_Z2_bars(results, H_grid, eta_grid, save_path=None):
    """
    Figure Z2 — Zumbach scalar Z = mean(Z(5),Z(10),Z(20)) per model.
    Left panel: Z vs H (fixed eta = eta_grid[1]).
    Right panel: Z vs eta (fixed H = H_grid[1]).
    """
    eta_fixed=eta_grid[len(eta_grid)//2]; H_fixed=H_grid[len(H_grid)//2]
    marks=["o","s","^","D"]; lss=["-","--","-.",":"]

    fig,axes=plt.subplots(1,2,figsize=(12,5))
    fig.suptitle(
        "Zumbach Scalar $Z = \\frac{1}{3}[Z(5)+Z(10)+Z(20)]$ "
        "— PDV-GL DGP (ground truth)",fontsize=11,fontweight="bold")

    # Left: Z vs H
    ax=axes[0]
    for ni,name in enumerate(ALL_NAMES):
        vals=[float(np.nanmean([s["zumbach"] for s in results[(H,eta_fixed)][name]]))
              if (H,eta_fixed) in results else np.nan for H in H_grid]
        ax.plot(H_grid,vals,color=_col(name),lw=2.,marker=marks[ni],
                ms=8,ls=lss[ni],label=name)
    ax.axhline(0,color=GRAY,lw=0.8,ls="--")
    ax.set_xticks(H_grid)
    ax.set_xlabel("Hurst $H$",fontsize=10)
    ax.set_ylabel("Zumbach scalar $Z$",fontsize=10)
    ax.set_title(f"$Z$ (single point: Table 8, Parameter Set 1)",fontsize=10)
    ax.legend(fontsize=9); ax.grid(True,alpha=0.22)
    for H in H_grid:
        if (H,eta_fixed) in results:
            dgpZ=float(np.nanmean([s["zumbach"] for s in results[(H,eta_fixed)]["DGP (PDV-GL)"]]))
            ax.annotate(f"{dgpZ:.3f}",(H,dgpZ),textcoords="offset points",
                        xytext=(0,6),ha="center",fontsize=7,color=GRAY)

    # Right: Z vs eta
    ax=axes[1]
    for ni,name in enumerate(ALL_NAMES):
        vals=[float(np.nanmean([s["zumbach"] for s in results[(H_fixed,eta)][name]]))
              if (H_fixed,eta) in results else np.nan for eta in eta_grid]
        ax.plot(eta_grid,vals,color=_col(name),lw=2.,marker=marks[ni],
                ms=8,ls=lss[ni],label=name)
    ax.axhline(0,color=GRAY,lw=0.8,ls="--")
    ax.set_xticks(eta_grid)
    ax.set_xlabel("Weight scale $w$",fontsize=10)
    ax.set_ylabel("Zumbach scalar $Z$",fontsize=10)
    ax.set_title(f"$Z$ vs $w$  (fixed $\\alpha={H_fixed}$)",fontsize=10)
    ax.legend(fontsize=9); ax.grid(True,alpha=0.22)

    plt.tight_layout(rect=[0,0,1,0.94])
    if save_path:
        fig.savefig(save_path,dpi=150,bbox_inches="tight")
        print(f"  Saved -> {save_path}")
    return fig


# ── Figure Z3: DGP diagnostics (returns + leverage + vol ACF) ─

def plot_Z3_diagnostics(results, H_grid, eta_fixed, save_path=None):
    """
    Figure Z3 — DGP realism diagnostics.
    One column per H. Four rows:
      0. Return histogram (log-y) with Gaussian overlay.
      1. Hill tail-index curve.
      2. Leverage ACF: Corr(r_t, V_{t+L}) vs L.
      3. Volatility ACF: Corr(V_t, V_{t+L}) vs L.
    """
    Hcols=[H for H in H_grid if (H,eta_fixed) in results]
    n_c=len(Hcols); lags=np.arange(1,21)

    fig,axes=plt.subplots(4,n_c,figsize=(5*n_c,14))
    if n_c==1: axes=axes[:,None]

    fig.suptitle(
        f"DGP + Model Diagnostics — PDV-GL (ground truth)\n"
        "Histogram | Hill | Leverage ACF | Volatility ACF",
        fontsize=11,fontweight="bold")

    for ci,H in enumerate(Hcols):
        data=results[(H,eta_fixed)]
        axes[0,ci].set_title(f"$\\alpha={H}$",fontsize=11,fontweight="bold")

        for name in ALL_NAMES:
            r_all=np.concatenate([s["logret"] for s in data[name]])
            col=_col(name); lw_=2.2 if "DGP" in name else 1.5

            # Row 0: histogram
            ax=axes[0,ci]
            lo,hi=np.percentile(r_all,0.3),np.percentile(r_all,99.7)
            ax.hist(r_all,bins=100,density=True,histtype="step",lw=lw_,color=col,label=name)
            if "DGP" in name:
                xr=np.linspace(lo,hi,300)
                ax.plot(xr,stats.norm.pdf(xr,0,r_all.std()),
                        color=GRAY,lw=0.9,ls=":",label="Gaussian ref")
            ax.set_yscale("log"); ax.set_xlim(lo,hi)
            if ci==0: ax.set_ylabel("Density (log)",fontsize=9)
            ax.legend(fontsize=6); ax.grid(True,alpha=0.18)

            # Row 1: Hill
            ax=axes[1,ci]
            ks,al=hill_curve(r_all)
            ax.plot(ks,al,color=col,lw=lw_,label=name)
            if ci==0: ax.set_ylabel(r"Hill $\hat\alpha$",fontsize=9)
            ax.legend(fontsize=6); ax.grid(True,alpha=0.18)

            # Row 2: Leverage ACF
            ax=axes[2,ci]
            lev_m=np.nanmean([s["lev_vec"] for s in data[name]],axis=0)
            ax.plot(lags,lev_m,color=col,lw=lw_,label=name)
            if ci==0: ax.set_ylabel("Leverage ACF\n$\\mathrm{Corr}(r_t,V_{t+L})$",fontsize=8)
            ax.axhline(0,color=GRAY,lw=0.6,ls="--")
            ax.legend(fontsize=6); ax.grid(True,alpha=0.18)

            # Row 3: Vol ACF
            ax=axes[3,ci]
            vacf_m=np.nanmean([s["vol_acf_vec"] for s in data[name]],axis=0)
            ax.plot(lags,vacf_m,color=col,lw=lw_,label=name)
            if ci==0: ax.set_ylabel("Vol ACF\n$\\mathrm{Corr}(V_t,V_{t+L})$",fontsize=8)
            ax.set_xlabel("Lag $L$ (days)",fontsize=9)
            ax.axhline(0,color=GRAY,lw=0.6,ls="--")
            ax.legend(fontsize=6); ax.grid(True,alpha=0.18)

    plt.tight_layout(rect=[0,0,1,0.94])
    if save_path:
        fig.savefig(save_path,dpi=150,bbox_inches="tight")
        print(f"  Saved -> {save_path}")
    return fig


# ── Figure Z4: Summary across H and eta grids ────────────────

def plot_Z4_summary(results, H_grid, eta_grid, save_path=None):
    """
    Figure Z4 — Full summary grid (4 rows).
      0. Zumbach scalar Z vs H  (cols = eta values)
      1. Leverage mean vs H
      2. H_hat vs H
      3. Score S vs H
    """
    marks=["o","s","^","D"]; lss=["-","--","-.",":"]
    cols_eta=eta_grid
    n_eta=len(cols_eta)

    fig,axes=plt.subplots(4,n_eta,figsize=(5*n_eta,14),sharey="row")
    if n_eta==1: axes=axes[:,None]
    fig.suptitle(
        "Summary — PDV-GL DGP (ground truth)"
        r"  |  $V_t=\xi_0\exp(\eta W_t^H-\eta^2t^{2H}/2)$",
        fontsize=11,fontweight="bold")

    for ci,eta in enumerate(cols_eta):
        Hok=[H for H in H_grid if (H,eta) in results]
        if not Hok: continue

        axes[0,ci].set_title("",fontsize=11,fontweight="bold")

        for ni,name in enumerate(ALL_NAMES):
            col=_col(name); mk=marks[ni]; ls=lss[ni]

            # Row 0: Zumbach
            ax=axes[0,ci]
            Z_=[float(np.nanmean([s["zumbach"] for s in results[(H,eta)][name]])) for H in Hok]
            ax.plot(Hok,Z_,color=col,lw=1.8,marker=mk,ms=6,ls=ls,label=name)
            ax.axhline(0,color=GRAY,lw=0.6,ls="--")
            if ci==0: ax.set_ylabel("Zumbach $Z$",fontsize=10)
            ax.set_xticks(Hok); ax.legend(fontsize=7); ax.grid(True,alpha=0.20)

            # Row 1: Leverage
            ax=axes[1,ci]
            L_=[float(np.nanmean([s["leverage"] for s in results[(H,eta)][name]])) for H in Hok]
            ax.plot(Hok,L_,color=col,lw=1.8,marker=mk,ms=6,ls=ls,label=name)
            ax.axhline(0,color=GRAY,lw=0.6,ls="--")
            if ci==0: ax.set_ylabel("Leverage $L$",fontsize=10)
            ax.set_xticks(Hok); ax.legend(fontsize=7); ax.grid(True,alpha=0.20)

            # Row 2: H_hat
            ax=axes[2,ci]
            H_=[float(np.nanmean([s["H_hat"] for s in results[(H,eta)][name]])) for H in Hok]
            ax.plot(Hok,H_,color=col,lw=1.8,marker=mk,ms=6,ls=ls,label=name)
        ax=axes[2,ci]
        ax.plot(Hok,Hok,color=CORAL,lw=1.2,ls=":",label="Identity")
        if ci==0: ax.set_ylabel(r"$\hat H$",fontsize=10)
        ax.set_xticks(Hok); ax.legend(fontsize=7); ax.grid(True,alpha=0.20)

        # Row 3: Score
        ax=axes[3,ci]
        for ni,name in enumerate(ALL_NAMES):
            col=_col(name); mk=marks[ni]; ls=lss[ni]
            S_=[float(np.nanmean([score_fn(s,results[(H,eta)]["ref"])
                                  for s in results[(H,eta)][name]])) for H in Hok]
            ax.plot(Hok,S_,color=col,lw=1.8,marker=mk,ms=6,ls=ls,label=name)
        ax.axhline(0,color=GRAY,lw=0.6,ls="--")
        if ci==0: ax.set_ylabel("Score $S$",fontsize=10)
        ax.set_xlabel("Hurst $H$",fontsize=10)
        ax.set_xticks(Hok); ax.legend(fontsize=7); ax.grid(True,alpha=0.20)

    row_titles=["Zumbach $Z$","Leverage $L$",r"$\hat H$","Score $S$"]
    for ri,tt in enumerate(row_titles):
        axes[ri,0].set_title(tt,fontsize=9,loc="left",pad=2)

    plt.tight_layout(rect=[0,0,1,0.95])
    if save_path:
        fig.savefig(save_path,dpi=150,bbox_inches="tight")
        print(f"  Saved -> {save_path}")
    return fig


def plot_leverage_figure(results, H_grid, eta_fixed,
                         sel_lags=(1, 5, 10, 20), save_path=None):
    """
    Figure L1 — Leverage effect: Corr(r_t, V_{t+L}) vs lag L.
    Reuses the per-path lev_vec already computed in compute_statistics
    (lags 1..20 days). Conventionally NEGATIVE in equity markets: a
    negative return today raises future variance (the classical
    "leverage effect", distinct from -- but related in spirit to --
    the Zumbach time-reversal asymmetry of the other figures).

    One column per H value (fixed eta). Three rows:
      0. Leverage curve: Gaussian... mean ± 1std band, per source.
      1. Leverage at selected lags (grouped bar chart).
      2. Leverage(1d) ranked horizontal bar (most negative = strongest
         leverage effect, at the bottom).
    """
    lags = np.arange(1, 21)
    Hcols = [H for H in H_grid if (H, eta_fixed) in results]
    n_c = len(Hcols)

    fig, axes = plt.subplots(3, n_c, figsize=(6*n_c, 13.5))
    if n_c == 1: axes = axes[:, None]

    fig.suptitle(
        "Leverage Effect — PDV-GL DGP (ground truth)\n"
        r"$L(L_{\rm lag}) = \mathrm{Corr}(r_t,\,V_{t+L_{\rm lag}})$"
        r"  —  $\mathbf{negative}$ = classical leverage effect",
        fontsize=12, fontweight="bold")

    for ci, H in enumerate(Hcols):
        data = results[(H, eta_fixed)]

        curves = {}
        for name in ALL_NAMES:
            arr = np.array([s["lev_vec"] for s in data[name]], dtype=float)
            curves[name] = dict(m=np.nanmean(arr, axis=0),
                                s=np.nanstd(arr, axis=0))

        col_title = "PDV-GL Table 8, Parameter Set 1"

        # ── Row 0: leverage curve ──────────────────────────────────
        ax = axes[0, ci]
        for name in ALL_NAMES:
            c = curves[name]; col = _col(name)
            lw = 2.4 if "DGP" in name else 1.6
            ax.fill_between(lags, c["m"]-c["s"], c["m"]+c["s"],
                            alpha=0.13, color=col)
            ax.plot(lags, c["m"], color=col, lw=lw, marker="o", ms=3,
                    label=name, zorder=4)
        ax.axhline(0, color=GRAY, lw=0.9, ls="--", alpha=0.6)
        ax.set_xlabel("Lag (days)", fontsize=9)
        if ci == 0: ax.set_ylabel(r"$\mathrm{Corr}(r_t,V_{t+L})$", fontsize=10)
        ax.set_title(f"{col_title}\nLeverage curve\n(most negative at short lag)",
                     fontsize=9)
        ax.legend(fontsize=7, loc="lower right")
        ax.grid(True, alpha=0.18)

        # ── Row 1: leverage at selected lags (grouped bar) ─────────
        ax = axes[1, ci]
        n_src = len(ALL_NAMES); bw = 0.70 / n_src
        x_pos = np.arange(len(sel_lags))
        for si, name in enumerate(ALL_NAMES):
            c = curves[name]; col = _col(name)
            vals = [float(c["m"][L-1]) for L in sel_lags]
            offset = (si - n_src/2 + 0.5) * bw
            ax.bar(x_pos + offset, vals, bw*0.92,
                   color=col, alpha=0.85, label=name)
        ax.axhline(0, color=GRAY, lw=0.9, ls="--")
        ax.set_xticks(x_pos)
        ax.set_xticklabels([f"{L}d" for L in sel_lags])
        if ci == 0: ax.set_ylabel("Leverage (mean)", fontsize=10)
        ax.set_title("Leverage at selected lags", fontsize=9)
        ax.legend(fontsize=6); ax.grid(True, axis="y", alpha=0.18)

        # ── Row 2: Leverage(1d) ranked horizontal bar ──────────────
        ax = axes[2, ci]
        l1d = {name: float(curves[name]["m"][0]) for name in ALL_NAMES}
        sorted_items = sorted(l1d.items(), key=lambda x: x[1])  # most negative first
        names_s = [it[0] for it in sorted_items]
        vals_s  = [it[1] for it in sorted_items]
        cols_s  = [_col(n) for n in names_s]
        y_ = np.arange(len(names_s))
        ax.barh(y_, vals_s, color=cols_s, alpha=0.85)
        ax.axvline(0, color=GRAY, lw=0.8, ls="--")
        ax.set_yticks(y_); ax.set_yticklabels(names_s, fontsize=8)
        if ci == 0: ax.set_xlabel("Leverage at lag = 1 day", fontsize=9)
        ax.set_title("Leverage(1d) ranked\nmore negative = stronger effect",
                     fontsize=8)
        ax.grid(True, axis="x", alpha=0.18)
        for yi, (name, v) in enumerate(zip(names_s, vals_s)):
            ha = "right" if v < 0 else "left"
            pad = -0.004 if v < 0 else 0.004
            ax.text(v + pad, yi, f"{v:.3f}", ha=ha, va="center", fontsize=7)

    plt.tight_layout(rect=[0, 0, 1, 0.92])
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
{\large Zumbach Effect: ESN vs.\ the PDV-GL Model Used as Ground-Truth DGP}\\[0.2cm]
{\normalsize Othmane Zarhali --- Paris Dauphine / CNRS}
\end{center}
\vspace{0.15cm}\noindent\rule{\textwidth}{0.8pt}

\begin{abstract}
We evaluate the Zumbach time-reversal asymmetry effect produced by the
ESN A2-981003 and Quadratic Rough Heston (QRH) models when the
\textbf{ground-truth data-generating process (DGP) is itself the
Guyon \& Lekeufack path-dependent volatility model (PDV-GL)}, in its
\textbf{4-factor Markovian} form \cite{guyon2023}. The DGP uses the
literature-calibrated ``Table 8, Parameter Set 1 (Realistic paths)''
values ($\beta_0=0.04,\beta_1=-0.13,\beta_2=0.65,\beta_{1,2}=0$,
$\lambda_{1,0}=55,\lambda_{1,1}=10,\theta_1=0.25$,
$\lambda_{2,0}=20,\lambda_{2,1}=3,\theta_2=0.5$), rescaled post-hoc to
the common $20\%$ annualised target (\S1) -- a single, fully-specified
calibration, with no swept grid parameter.
\textbf{Only the ESN and QRH are calibrated and compared against this
DGP} --- a separately-calibrated copy of PDV-GL itself is deliberately
not run, since PDV-GL is now the ground truth and comparing it against
a freshly-fit copy of itself would add little. The ESN is run at its
current optimal architecture ($N_r{=}96,N_z{=}12$); both models are
calibrated using the full 11-term data-adaptive score $S$
(Eq.~\eqref{eq:score}) with the \textbf{smooth Gaussian-kernel}
proximity. The ESN's read-out combines a linear rough factor, a
per-mode $z$-readout profile, and a symmetric quadratic-term matrix,
all calibrated per DGP scenario alongside $H$ and the reservoir rate
bounds (\S\ref{sec:calib}). The Zumbach effect is assessed via the
per-lag curve $Z(L)$ and the scalar
$Z=\frac{1}{3}[Z(5)+Z(10)+Z(20)]$, supplemented by leverage ACF, vol ACF,
and return tail diagnostics.
\end{abstract}

\tableofcontents\vspace{0.4cm}

\section{DGP: PDV-GL as Ground Truth}

The data-generating process is the Guyon \& Lekeufack path-dependent
volatility model \cite{guyon2023}, in its \textbf{4-factor Markovian}
form -- the standard finite-dimensional approximation obtained by
replacing each power-law memory kernel with a mixture of two
exponential (OU-type) kernels, making the model fast to simulate as a
genuine Markov process:
\begin{equation}
  \frac{dS_t}{S_t} = \sigma_t\,dW_t,
  \label{eq:pdvgl_ret}
\end{equation}
\begin{equation}
  \sigma_t = \beta_0 + \beta_1 R_{1,t} + \beta_2\sqrt{R_{2,t}} + \beta_{1,2}\,R_{1,t}\sqrt{R_{2,t}},
  \label{eq:pdvgl_var}
\end{equation}
\begin{equation}
  R_{1,t} = \theta_1 X_{1,0,t} + (1-\theta_1) X_{1,1,t},
  \qquad
  R_{2,t} = \theta_2 X_{2,0,t} + (1-\theta_2) X_{2,1,t},
  \label{eq:pdvgl_factors}
\end{equation}
where $X_{1,j,t}$ is an exponentially-weighted moving average of past
\emph{signed} returns with mean-reversion rate $\lambda_{1,j}$ (the
\textbf{leverage} factor, a mixture of a fast and a slow exponential
kernel), and $X_{2,j,t}$ is an EWMA of past \emph{squared} returns with
rate $\lambda_{2,j}$ (the \textbf{magnitude} factor). Unlike every
other DGP in this line of protocols, PDV-GL is used here as
\textbf{ground truth}, not as a candidate model.

\paragraph{"True" parameters -- literature-calibrated Table 8,
Parameter Set 1 (Realistic paths).}
\begin{table}[h]
\centering
\caption{PDV-GL 4-factor Markovian ground-truth parameters (Table 8, Parameter Set 1).}
\label{tab:pdvgl4f}
\begin{tabular}{@{}lr@{}}
\toprule
Parameter & Value \\
\midrule
$\beta_0$ & 0.04 \\
$\beta_1$ & $-0.13$ \\
$\beta_2$ & 0.65 \\
$\beta_{1,2}$ & 0 (unused for Set 1) \\
$\lambda_{1,0}$ & 55 \\
$\lambda_{1,1}$ & 10 \\
$\theta_1$ & 0.25 \\
$\lambda_{2,0}$ & 20 \\
$\lambda_{2,1}$ & 3 \\
$\theta_2$ & 0.5 \\
\bottomrule
\end{tabular}
\end{table}
$\lambda$ rates are in annualised units; each is applied at the daily
step via $e^{-\lambda\Delta t}$, $\Delta t=1/252$.

\paragraph{Rescaling.}
Run literally with these parameters (and daily returns in decimal
units), the model's realised annualised volatility comes out at only
$\approx4\%$ -- the source table's implicit return-scale convention
could not be verified from the parameters alone. A single, documented
linear rescaling $X\mapsto X\cdot(\mathrm{TV\_DAY}/\widehat\sigma)$ is
applied post-hoc to bring the level to the same $20\%$ annualised target
used throughout this line of protocols. This is a pure amplitude change:
every \emph{normalised} statistic used in this protocol ($\hat H$,
kurtosis, leverage correlation, and the $Z(\tau)/\sigma^3$ cross-moment
ratio) is left \textbf{exactly unchanged}, since numerator and
denominator scale identically under a uniform linear rescaling.

\paragraph{Zumbach effect mechanism.}
Because $\beta_1<0$, a large \emph{negative} past return raises
$R_{1,t}$ with the sign needed to raise $\sigma_t$ and hence future
variance -- a genuine, signed leverage asymmetry. At these literature
parameters, the realised leverage correlation is $\approx-0.04$ to
$-0.17$ depending on the exact path sample, and the windowed Zumbach
scalar $Z$ (Eq.~\ref{eq:score}'s calibration statistic) comes out
strongly positive ($Z\approx0.22$ in the run reported here) -- a
clearly non-trivial, sustained effect rather than a borderline one.

\paragraph{Simulation.}
Eq.~\eqref{eq:pdvgl_factors} is iterated as an exact discretisation of
the four OU-type factors ($X_{1,0},X_{1,1},X_{2,0},X_{2,1}$) at the
daily step, using $e^{-\lambda_{i,j}\Delta t}$ decay -- a genuine Markov
process in $(X_{1,0},X_{1,1},X_{2,0},X_{2,1})$, no path history storage
needed (unlike the power-law-kernel DGPs used elsewhere in this line of
protocols).

\paragraph{Parameter grid.}
None -- Parameter Set 1 is a single, fully-specified literature
calibration (Table~\ref{tab:pdvgl4f}); every path uses the same fixed
parameters.

\section{Models}
\subsection{ESN A2-981003 (current optimal architecture)}

The ESN uses Gaussian Brownian innovations, with 7 shared
hyperparameters and $N_r=96,N_z=12$. \texttt{rough\_scale} and the
$z$-readout ($\texttt{zr\_lo},\texttt{zr\_hi}$, plus the symmetric
quadratic matrix coefficients $(m_1,m_2)$) are calibrated per DGP cell
(\S\ref{sec:calib}), not part of this fixed table.
\begin{table}[h]
\centering
\caption{The 7 shared hyperparameters of the optimal architecture.}
\label{tab:optarch}
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
\midrule
\texttt{rough\_scale} & \textbf{calibrated per DGP cell} (\S\ref{sec:calib}) \\
\texttt{zr\_lo}, \texttt{zr\_hi} & \textbf{calibrated per DGP cell} (was a single fixed \texttt{z\_readout}=0.0516) \\
$\kappa_{\rm quad}$ & \textbf{calibrated per DGP cell} (new) \\
\bottomrule
\end{tabular}
\end{table}
Single $\varepsilon_t\sim\mathcal{N}(0,1)$ drives both the Ito log-return
and the diagonal OU rough bank. The negative instantaneous leverage
$\kappa_0=-0.434991<0$ (admissibility, guaranteed by construction) gives
the ESN an endogenous Zumbach-like asymmetry -- structurally different
from PDV-GL's own path-dependent-feedback mechanism (\S1), so a good
match here is evidence the ESN's mechanism can \emph{mimic} the effect,
not that it replicates the same causal channel. $H$, the $r$-layer rate
window $[\lambda_{\rm lo},\lambda_{\rm hi}]$, and the $z$-layer rate
window $[a_{\rm lo},a_{\rm hi}]$ are \textbf{calibrated per $\alpha_r$
DGP cell} (\S\ref{sec:calib}) -- they are not part of the shared
architecture above.

\subsection{Quadratic Rough Heston (QRH)}
Bourgey \& Gatheral: $V_t=Y_t^2+c$, gamma kernel
$\kappa(\tau)=\nu/\Gamma(\alpha)\tau^{\alpha-1}e^{-\lambda\tau}$,
$\alpha=H+\tfrac{1}{2}$, Euler-Volterra with 252-step ring buffer.
Included alongside the ESN as a second, structurally different
candidate model (a genuine rough-volatility SDE, vs.\ the ESN's
Markovian reservoir approximation) against the same PDV-GL ground truth.

\section{Calibration: Full 11-Term Score, Smooth Kernel}
\label{sec:calib}

The data-adaptive score (identical structure to every other protocol in
this series, Eq.~\eqref{eq:score}), with the \textbf{smooth
Gaussian-kernel} proximity (protocol V5) rather than the original
piecewise-linear tent:
\begin{equation}
\begin{aligned}
S &= 2.2\,\bigl(1-f(\hat H)\bigr)+1.4\,\bigl(1-f(\bar\sigma)\bigr)+0.6\,\bigl(1-g(q_{995})\bigr)+0.2\,\bigl(1-g(\sigma_{\max})\bigr)\\
  &+0.8\,\bigl(1-f(\bar V_{\mathrm{ACF}})\bigr)+1.0\,\bigl(1-f(G_T)\bigr)+0.8\,\bigl(1-f(F_T)\bigr)+\mathbf{3.5}\,\bigl(1-f(Z)\bigr)\\
  &+1.1\,\bigl(1-f(L)\bigr)+0.7\,\bigl(1-f(K)\bigr)+0.5\,\bigl(1-f(A)\bigr)+5\cdot\mathrm{Stress},
\end{aligned}
\label{eq:score}
\end{equation}
where $f(x,c,s)=\exp(-\tfrac12((x-c)/s)^2)$,
$g(x,c,s)=\exp(-\tfrac12((x-c)_+/s)^2)$, both in $(0,1]$. $S\geq0$ is
\textbf{minimised directly} (no sign flip); $S=0$ is a perfect match on
all 11 stylised facts and no stress event; \textbf{$S=17.8$ is
worst-possible} (raised from $15.3$ by the Zumbach-weight increase
below). All centres $c_k$ and tolerances $s_k$ are DGP cross-path mean
and $\max(\mathrm{std},0.30|c_k|,\mathrm{floor}_k)$.

\paragraph{Zumbach score weighting.} The Zumbach weight is raised
from $1.0$ to $\mathbf{3.5}$ -- now the single largest weight in $S$,
exceeding even $\hat H$'s $2.2$ -- and QRH's calibration uses a
2-point multi-start. \S\ref{sec:zumbach-results} reports the full
root-cause diagnosis and result.

\paragraph{ESN read-out.} The ESN's read-out combines a linear rough
factor, a per-mode $z$-readout profile
$(\texttt{zr\_lo},\texttt{zr\_hi})$, and a symmetric quadratic-term
matrix $Q=(1/N_r)(m_1I+m_2qq^\top)$ added to the pre-activation. The
ESN calibrates its full inner-loop vector
\[
H,\ \lambda_{\rm lo},\ \lambda_{\rm hi},\ a_{\rm lo},\ a_{\rm hi},\
\Delta b_0,\ \mathrm{scale},\ \texttt{rough\_scale},\
\texttt{zr\_lo},\ \texttt{zr\_hi},\ m_1,\ m_2
\]
-- \textbf{12 dimensions} -- via an 11-point multi-start
Nelder--Mead scheme (rough/fast, persistent/slow, intermediate,
wide-band/multiscale, long-memory extreme, Zumbach
term1/term2-balanced, strong-amplitude, weak-amplitude,
isotropic-positive, isotropic-negative, and $q$-aligned-positive
starts), keeping the best-of-eleven result; only the 7 shared
architecture hyperparameters of Table~\ref{tab:optarch} stay fixed.

Calibrated parameters:
\begin{center}\renewcommand{\arraystretch}{1.3}
\begin{tabular}{@{}p{3.0cm}p{6.4cm}p{5.0cm}@{}}\toprule
Model & Parameters calibrated & Remark \\\midrule
ESN A2-981003 & $H$, $\lambda_{\rm lo}$, $\lambda_{\rm hi}$, $a_{\rm lo}$,
  $a_{\rm hi}$, $\Delta b_0$, scale $s$, \texttt{rough\_scale},
  \texttt{zr\_lo}, \texttt{zr\_hi}, $m_1$, $m_2$ &
  7 shared hyperparameters fixed; full 12-dim.\ inner-loop calibrated
  per $\alpha_r$ cell (11-start Nelder--Mead) \\
QRH & $H,\nu_{\rm vol},\lambda,c_{\rm frac}$ & 2-start Nelder--Mead
  (added a low-$H$ start; see \S\ref{sec:zumbach-results}) \\
\bottomrule\end{tabular}\end{center}

\section{Why Didn't QRH's Zumbach Fit, and the Fix}
\label{sec:zumbach-results}

\begin{table}[h]
\centering
\caption{Zumbach scalar $Z$, before vs.\ after the fix. DGP target: $Z=0.220$.}
\begin{tabular}{@{}lcc@{}}
\toprule
Model & $Z$, before & $Z$, after \\
\midrule
ESN A2-981003 & 0.173 & 0.182 \\
QRH & $-0.044$ (wrong sign) & \textbf{0.230} \\
\bottomrule
\end{tabular}
\end{table}

\paragraph{The problem.} Without the fix below, QRH's calibrated Zumbach
scalar comes out \emph{negative} ($-0.044$) against
the DGP's clearly positive $0.220$ -- not merely a weak match, but the
wrong sign entirely, plainly visible as the large grey-vs-blue gap in
the main Zumbach figure's difference column.

\paragraph{Root-cause diagnosis.} A dedicated sensitivity sweep of
QRH's four free parameters $(H,\nu_{\rm vol},\lambda,c_{\rm frac})$
revealed that QRH's Zumbach-generating capacity is strongly
\textbf{anti-correlated with its roughness $H$}: at low $H$
($\approx0.03$--$0.10$), QRH produces a robust, DGP-matching positive
Zumbach ($0.14$--$0.24$ across 15-path averages); at higher $H$
($\approx0.4$, needed to match the DGP's own $\hat H_c\approx0.32$),
QRH's Zumbach is weak and \textbf{noise-dominated} -- a 20-path test at
the actually-calibrated $H=0.42$ gave individual-path values ranging
from $-0.27$ to $+0.08$ (mean $-0.02$, std $0.085$: statistically
indistinguishable from zero, let alone matching the DGP).

The calibration was choosing $H\approx0.42$ because the score's $\hat H$
term (weight $2.2$) pulls toward matching the DGP's own $\hat
H_c\approx0.32$, while the (old) Zumbach weight of only $1.0$ was not
strong enough to outweigh that pull once QRH's Zumbach was already deep
in its noise floor. A back-of-envelope check confirms the sign of the
trade-off: moving from $H\approx0.42$ to $H\approx0.10$ costs roughly
$+1.1$ in the $\hat H$ term (further from target) but saves
$\approx1.0$ in the Zumbach term \emph{at the old weight of 1.0} --
close to a wash, so Nelder--Mead had little incentive to move; at the
new weight of $3.5$ the same Zumbach saving is worth $\approx3.0$,
clearly outweighing the $\hat H$ cost.

\paragraph{The fix, in two parts.} (1) Raise the Zumbach weight from
$1.0$ to $3.5$ (\S\ref{sec:calib}) so the trade-off above favours
matching Zumbach. (2) This alone was not sufficient: with QRH's single
calibration start biased toward the DGP's $\hat H_c$, Nelder--Mead
still converged to $H\approx0.41$ and $Z=-0.028$ -- barely moved. QRH's
calibration was therefore widened to a \textbf{2-point multi-start},
adding a dedicated low-$H$ starting point ($H_0=0.06$) alongside the
original $\hat H$-biased one, keeping the best-of-two result.

\paragraph{Result.} With both changes, QRH's calibration converges to
$H=0.09$, giving $Z=0.230$ -- matching the DGP's $0.220$ almost exactly,
a complete reversal from the wrong-signed $-0.044$ before. This confirms
the mismatch was a genuine calibration failure (a score trade-off
compounded by a biased single starting point), not a structural
limitation of QRH itself -- unlike the PDV-GL findings in the
Taylor-effect protocol, where no recalibration could produce the
missing effect. The ESN's Zumbach fit was already reasonable before this
fix ($0.173$) and improves slightly to $0.182$, since only the shared
score changed (its own architecture and calibration were not modified
for this fix).

\subsection{A second, deeper problem: matching the difference is not
enough}
\label{sec:zumbach-terms}

Even after the fix above, closer inspection of the \textbf{individual}
terms $\mathrm{Term\,1}=\mathrm{Corr}(\mathrm{past}\,R^2,\mathrm{fut}\,V)$
and $\mathrm{Term\,2}=\mathrm{Corr}(\mathrm{past}\,V,\mathrm{fut}\,R^2)$
(the two left columns of the main Zumbach figure) revealed that both
ESN and QRH matched the DGP's \emph{difference} $Z$ well while getting
\textbf{both individual terms wrong in the same direction}: the DGP has
$\mathrm{Term\,1}_c\approx0.221$ and a near-zero
$\mathrm{Term\,2}_c\approx0.011$ (i.e.\ essentially no
past-vol-predicts-future-magnitude channel at all), whereas both models
showed $\mathrm{Term\,1}\approx0.30$--$0.45$ \emph{and}
$\mathrm{Term\,2}\approx0.16$--$0.33$ -- both terms inflated together,
so their difference still landed near the DGP's $Z$ by cancellation,
for the wrong reason.

\begin{table}[h]
\centering
\caption{Individual terms, before vs.\ after this second fix. DGP target: Term1$=0.221$, Term2$=0.011$.}
\begin{tabular}{@{}lcccc@{}}
\toprule
 & \multicolumn{2}{c}{Term 1} & \multicolumn{2}{c}{Term 2}\\
Model & before & after & before & after \\
\midrule
ESN A2-981003 & 0.333 & 0.134 & 0.163 & $-0.007$ \\
QRH & 0.348 & 0.264 & 0.331 & \phantom{$-$}0.002 \\
\bottomrule
\end{tabular}
\end{table}

\paragraph{The fix.} The single Zumbach-difference score term
(\S\ref{sec:calib}, weight $3.5$) is replaced by \textbf{two} terms
matching $\mathrm{Term\,1}$ and $\mathrm{Term\,2}$
\emph{individually} (weight $1.75$ each, same combined weight, same
worst-possible $S=17.8$) -- a model can no longer satisfy the score by
getting the difference right while both channels are wrong. This alone
changed the optimisation landscape enough that QRH's own multi-start
needed a third, dedicated starting point (low $H$, low
$\nu_{\rm vol}=0.15$) to reliably find a configuration balancing both
terms simultaneously, and the ESN's multi-start needed a sixth start
(fast $r$-layer, slow-\emph{and}-strong $z$-layer:
$\lambda_{\rm hi}=5$, $a_{\rm lo},a_{\rm hi}=0.3,0.5$) for the same
reason -- both added alongside the existing starts, keeping the
overall best-of-set result.

\paragraph{Result.} $\mathrm{Term\,2}$ for both models drops from
$0.16$--$0.33$ to within $0.01$ of the DGP's near-zero target -- QRH's
own vol-of-vol $\nu_{\rm vol}$ collapses to $0.077$ (an order of
magnitude below its earlier calibrated value) to suppress this spurious
channel, and its kurtosis fit improves as a side effect ($6.8\to0.05$
excess kurtosis, much closer to the DGP's $0.008$). $\mathrm{Term\,1}$
for both models moves \emph{closer} to the DGP for QRH ($0.348\to0.264$)
but slightly \emph{past} it for the ESN ($0.333\to0.134$, now
undershooting rather than overshooting) -- a real, if partial,
remaining imperfection, visible as the blue curve sitting somewhat
below the grey one in the Term 1 panel. This is a smaller residual gap
than the original problem, and unlike the PDV-GL-as-DGP structural
limits found elsewhere in this line of protocols, both terms are now at
least of the right \emph{order of magnitude} and (for Term 2)
essentially exact for both models -- a substantially more faithful
Zumbach-effect reproduction than matching the difference alone provided.

\subsection{Individual-term fit at the current architecture}
\label{sec:zumbach-v8}

Matching only the aggregate difference $Z=\mathrm{Term1}-\mathrm{Term2}$
is not enough: a model can match $Z$ by having both terms elevated in
tandem. Table~\ref{tab:zumbachterms} reports the individual terms.

\begin{table}[h]
\centering
\caption{Individual Zumbach terms and aggregate $Z$. DGP target: Term1$=0.221$, Term2$=0.011$, $Z=0.220$.}
\label{tab:zumbachterms}
\begin{tabular}{@{}lccc@{}}
\toprule
Model & Term 1 & Term 2 & $Z$ \\
\midrule
ESN A2-981003 & 0.160 & 0.003 & 0.156 \\
QRH & 0.271 & 0.021 & 0.250 \\
\bottomrule
\end{tabular}
\end{table}

\textbf{Reading.} QRH matches both individual terms closely. The ESN
matches Term 2 well (close to the DGP's near-zero value) but
undershoots Term 1, giving an aggregate $Z$ that is directionally
correct but below target. This individual-term view is the reason the
Zumbach weight was raised and QRH's start widened (\S\ref{sec:calib});
matching the scalar $Z$ alone would not have revealed the Term 2
mismatch that a more casual calibration could otherwise paper over.

\section{Zumbach Diagnostic}
\label{sec:windowed}

\begin{equation}
  Z(L) = \mathrm{Corr}\!\bigl(\mathrm{past}_L R^2,\,\mathrm{fut}_L V\bigr)
        - \mathrm{Corr}\!\bigl(\mathrm{past}_L V,\,\mathrm{fut}_L R^2\bigr),
  \label{eq:zumbach_scalar}
\end{equation}
where $\mathrm{past}_L R^2=(X_t{-}X_{t-L})^2$ and
$\mathrm{fut}_L V=(V_{t+L}{-}V_t)/L$. $Z(L)>0$ means the
past-magnitude$\to$future-vol channel [$\mathrm{Corr}(\mathrm{past}_LR^2,\mathrm{fut}_LV)$]
dominates the reverse, time-reversed channel
[$\mathrm{Corr}(\mathrm{past}_LV,\mathrm{fut}_LR^2)$] -- the genuine
Zumbach / time-reversal-asymmetry signature. This is the SAME windowed,
correlation-based metric used for calibration (the scalar
$Z=\tfrac{1}{3}[Z(5){+}Z(10){+}Z(20)]$ enters the score with weight
1.0) and for Figure Z1. It depends only on magnitudes, not the sign of
returns, and so is well-defined even for DGPs with no leverage term at
all.

The ESN carries its own, structurally distinct source of asymmetry (the
fixed $\kappa_0<0$ leverage built into its architecture, \S2) and the
PDV-GL DGP its own path-dependent feedback (\S1) -- a good match
between the ESN's and the DGP's $Z(L)$ is therefore evidence that the
ESN's mechanism can \emph{mimic} the magnitude of the effect, not that
it reproduces the same causal channel.

\section{Figures}

\textbf{Main figure (Zumbach windowed metric).}
Three columns -- one per term of $Z(L)=\mathrm{Term\,1}-\mathrm{Term\,2}$
-- each overlaying all sources (DGP, ESN, QRH):
\begin{itemize}[itemsep=2pt]
  \item Column 1: $\mathrm{Corr}(\mathrm{past}\,R^2,\mathrm{fut}\,V)$
    (the genuine, forward Zumbach channel).
  \item Column 2: $\mathrm{Corr}(\mathrm{past}\,V,\mathrm{fut}\,R^2)$
    (the reverse, time-reversed channel).
  \item Column 3: $Z(L)=$ Column 1 $-$ Column 2 (the Zumbach metric itself).
\end{itemize}
Two rows per column: (i) curve vs lag, Gaussian-smoothed mean (line,
$\sigma_{\rm sm}=2$d) + dots + $\pm1$std band; (ii) grouped bars at
lags 1d, 5d, 10d, 20d. Showing the two terms side by side, rather than
overlaid in a single panel, makes it easier to see \emph{which} channel
drives any difference between sources -- e.g.\ whether a model matches
the DGP's $Z(L)$ because both terms individually match, or because two
mismatched terms happen to cancel to a similar difference.

\textbf{Figure Z1.} Windowed $Z(L)$ curves (row per source -- DGP, ESN, QRH -- column per $\alpha_r$).

\textbf{Figure Z2.} Scalar $Z$ vs $\alpha_r$, DGP vs ESN vs QRH.

\textbf{Figure Z3.} Return diagnostics: histogram | Hill | leverage ACF | vol ACF.

\textbf{Figure Z4.} Summary: $Z$, leverage, $\hat H$, score $S$ vs $\alpha_r$.

\textbf{Figure L1 (new).} Dedicated leverage-effect figure, one column
per $\alpha_r$, three rows:
\begin{enumerate}[itemsep=2pt]
  \item Leverage curve $\mathrm{Corr}(r_t,V_{t+L})$ vs lag $L$ (1--20
    days): mean $\pm1$std band, per source. Conventionally negative in
    equity markets -- distinct from, but related in spirit to, the
    Zumbach time-reversal asymmetry of the figures above.
  \item Grouped bars of the leverage curve at lags 1d, 5d, 10d, 20d.
  \item Horizontal bar of leverage at lag $=1$ day, ranked most negative
    at the bottom (strongest leverage effect); values annotated.
\end{enumerate}
Reuses the per-path \texttt{lev\_vec} array already computed in
\texttt{compute\_statistics} (the same statistic entering the score's
leverage term, weight 1.1) -- no new low-level computation, just a
dedicated, richer visualisation (mean$\pm$std band, selected-lag bars,
ranked comparison) than the single mean-line panel embedded in Figure Z3.

\section{Monte Carlo Setup}

\begin{center}\renewcommand{\arraystretch}{1.4}
\begin{tabular}{ll}\toprule Parameter & Value \\\midrule
DGP paths $n_{\rm dgp}$ & 30 \\
Calibration paths $n_{\rm cal}$ & 5 \\
Calibration days $T_{\rm cal}$ & 600 \\
Calibration burn-in & 150 \\
Final paths $n_{\rm sim}$ & 30 \\
Final days $T_{\rm sim}$ & 2500 \\
Final burn-in & 500 \\
$\Delta t$ & 1 (daily) \\
DGP & PDV-GL Table 8, Parameter Set 1 (fixed, no grid) \\
DGP $\beta_0,\beta_1,\beta_2$ & $0.04,\ -0.13,\ 0.65$ \\
DGP $\lambda_{1,0},\lambda_{1,1},\theta_1$ & $55,\ 10,\ 0.25$ \\
DGP $\lambda_{2,0},\lambda_{2,1},\theta_2$ & $20,\ 3,\ 0.5$ \\
Smoothing $\sigma$ & 2 days \\
\bottomrule\end{tabular}\end{center}

\begin{thebibliography}{9}
\bibitem{guyon2023} J.\ Guyon, J.\ Lekeufack.
\textit{Volatility is (mostly) path-dependent}. QF 23(9), 2023.
\bibitem{zumbach2010} G.\ Zumbach.
\textit{Time reversal invariance in finance}. Quantitative Finance
9(5):505--515, 2010.
\end{thebibliography}
\end{document}
"""

# ── Cross-moment Zumbach curve (exact figure format) ─────────

def zumbach_crossmoment_curve(st_list, max_lag=30):
    """
    Per-path cross-moment Zumbach curves.
    Z(tau)  = E[r^2_{t+tau} * r_t] - E[r_{t+tau} * r^2_t]
    Zpos(tau) = E[r^2_{t+tau} * r_t]   (future-squared × past)
    Zneg(tau) = E[r_{t+tau} * r^2_t]   (future × past-squared)

    Convention from the uploaded figure:
        Negative Z(tau) = Zumbach effect present
        (past squared returns predict HIGHER future vol more than vice versa,
         but the cross-moment Z = E[r^2_{t+tau} r_t] - E[r_{t+tau} r^2_t] < 0
         because with rho < 0 a large NEGATIVE r_t predicts high r^2_{t+tau})

    Returns:
        Z_m, Z_s    : cross-path mean and std of Z(tau),   shape (max_lag,)
        Zp_m, Zp_s  : cross-path mean and std of Zpos(tau)
        Zn_m, Zn_s  : cross-path mean and std of Zneg(tau)
    """
    lags = np.arange(1, max_lag + 1)
    Z_all=[]; Zp_all=[]; Zn_all=[]
    for st in st_list:
        r = st["logret"]
        r2 = r**2
        Zpos = np.array([np.mean(r2[k:] * r[:-k]) for k in lags])
        Zneg = np.array([np.mean(r[k:]  * r2[:-k]) for k in lags])
        Z_all.append(Zpos - Zneg)
        Zp_all.append(Zpos); Zn_all.append(Zneg)
    Z_a  = np.array(Z_all);  Zp_a = np.array(Zp_all); Zn_a = np.array(Zn_all)
    return (np.nanmean(Z_a, 0),  np.nanstd(Z_a, 0),
            np.nanmean(Zp_a, 0), np.nanstd(Zp_a, 0),
            np.nanmean(Zn_a, 0), np.nanstd(Zn_a, 0))


def plot_zumbach_figure(results, H_grid, eta_fixed,
                        max_lag=30, smooth_sigma=2.0,
                        sel_lags=(1, 5, 10, 20),
                        save_path=None):
    """
    THREE COLUMNS, one per term of the windowed Zumbach metric (replacing
    the earlier "one column per alpha_r" layout, now that there is only
    one alpha_r value to show):

        Z(L) = Corr(past_R^2, future_V) - Corr(past_V, future_R^2)
               \\_______ C1 _______/     \\_______ C2 _______/

    Column 0 — C1(L) = Corr(past_R^2, future_V)   (the genuine, forward
               Zumbach channel: past magnitude -> future vol)
    Column 1 — C2(L) = Corr(past_V, future_R^2)   (the reverse channel:
               past vol -> future magnitude)
    Column 2 — Z(L)  = C1(L) - C2(L)               (the Zumbach metric
               itself; POSITIVE = time-reversal-asymmetry present)

    Each column overlays ALL sources (DGP, ESN, QRH). Two rows:
      Row 0 — curve vs lag: Gaussian-smoothed mean + unsmoothed dots +
              ±1std shaded band.
      Row 1 — grouped bar chart at selected lags (1d, 5d, 10d, 20d).
    """
    from scipy.ndimage import gaussian_filter1d

    Hcols = [H for H in H_grid if (H, eta_fixed) in results]
    if not Hcols:
        return None
    H = Hcols[0]   # single fixed PDV-GL Table 8 Parameter Set 1 -- no
                   # alpha_r grid left to loop over (see module docstring)
    data = results[(H, eta_fixed)]
    lags = np.arange(1, max_lag + 1)

    # Pre-compute smoothed curves for all sources
    curves = {}
    for name in ALL_NAMES:
        Zm,Zs,C1m,C1s,C2m,C2s = zumbach_windowed_decomp(data[name], max_lag)
        curves[name] = dict(
            Zm=Zm, Zs=Zs, C1m=C1m, C1s=C1s, C2m=C2m, C2s=C2s,
            Zm_sm =gaussian_filter1d(np.nan_to_num(Zm),  smooth_sigma),
            C1m_sm=gaussian_filter1d(np.nan_to_num(C1m), smooth_sigma),
            C2m_sm=gaussian_filter1d(np.nan_to_num(C2m), smooth_sigma),
        )

    panels = [
        ("C1m", "C1s", "C1m_sm", r"$\mathrm{Corr}(\mathrm{past}\,R^2,\,\mathrm{fut}\,V)$",
         "Term 1: past $R^2$ $\\to$ future $V$", "upper right"),
        ("C2m", "C2s", "C2m_sm", r"$\mathrm{Corr}(\mathrm{past}\,V,\,\mathrm{fut}\,R^2)$",
         "Term 2: past $V$ $\\to$ future $R^2$", "upper right"),
        ("Zm",  "Zs",  "Zm_sm",  r"$Z(L) = \mathrm{Term\ 1} - \mathrm{Term\ 2}$",
         "Difference: Zumbach metric $Z(L)$", "upper right"),
    ]
    n_c = len(panels)

    fig, axes = plt.subplots(2, n_c, figsize=(6*n_c, 10))

    fig.suptitle(
        "Zumbach Time-Reversal Asymmetry — PDV-GL DGP (ground truth)\n"
        "PDV-GL Table 8, Parameter Set 1\n"
        r"$Z(L) = \mathrm{Corr}(\mathrm{past}\,R^2,\,\mathrm{fut}\,V) - \mathrm{Corr}(\mathrm{past}\,V,\,\mathrm{fut}\,R^2)$"
        r"  —  $\mathbf{positive}$ for all $L$ = Zumbach effect",
        fontsize=12, fontweight="bold")

    for pi, (mkey, skey, smkey, ylabel, title, loc) in enumerate(panels):
        # ── Row 0: curve vs lag ──────────────────────────────────────
        ax = axes[0, pi]
        for name in ALL_NAMES:
            c = curves[name]; col = _col(name)
            lw = 2.4 if "DGP" in name else 1.6
            ax.fill_between(lags, c[mkey]-c[skey], c[mkey]+c[skey],
                            alpha=0.13, color=col)
            ax.scatter(lags, c[mkey], color=col, s=7, alpha=0.35, zorder=3)
            ax.plot(lags, c[smkey], color=col, lw=lw, label=name, zorder=4)
        ax.axhline(0, color=GRAY, lw=0.9, ls="--", alpha=0.6)
        ax.set_xlabel("Lag $L$ (days)", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.legend(fontsize=7, loc=loc)
        ax.grid(True, alpha=0.18)

        # ── Row 1: bar chart at selected lags ─────────────────────────
        ax = axes[1, pi]
        n_src = len(ALL_NAMES); bw = 0.70 / n_src
        x_pos = np.arange(len(sel_lags))
        for si, name in enumerate(ALL_NAMES):
            c = curves[name]; col = _col(name)
            vals = [float(c[smkey][L-1]) for L in sel_lags]
            offset = (si - n_src/2 + 0.5) * bw
            ax.bar(x_pos + offset, vals, bw*0.92,
                   color=col, alpha=0.85, label=name)
        ax.axhline(0, color=GRAY, lw=0.9, ls="--")
        ax.set_xticks(x_pos)
        ax.set_xticklabels([f"{L}d" for L in sel_lags])
        ax.set_ylabel(f"{ylabel} (smoothed)", fontsize=9)
        ax.set_title("At selected lags", fontsize=9)
        ax.legend(fontsize=6); ax.grid(True, axis="y", alpha=0.18)

    plt.tight_layout(rect=[0, 0, 1, 0.90])
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"  Saved -> {save_path}")
    return fig


def write_latex(path):
    with open(path,"w") as f: f.write(LATEX.lstrip())
    print(f"  Saved -> {path}")

# ============================================================
# 8.  Main pipeline
# ============================================================

def run(H_grid=None, eta_grid=None, rho=RHO_DGP,
        n_dgp=30, n_sim=30, T=2500, burn=500,
        n_cal=5, T_cal=600, burn_cal=150, dt=1.0,
        eta_fixed=None, save_prefix="rough_bergomi_zumbach",
        latex_path=None):
    """
    NOTE: parameter names H_grid/eta_grid/eta_fixed/rho are kept for
    interface compatibility with earlier drafts of this protocol; here
    H_grid is actually the PDV-GL DGP's "true" short-memory/rough
    exponent grid (TRUE_ALPHA_R_GRID). eta_grid/eta_fixed/rho are unused
    (kept only so the CLI/run() signature doesn't need to change again).

    PDV-GL is now the GROUND-TRUTH DGP itself (not a calibrated model),
    and ONLY the ESN is calibrated against it and compared -- QRH and
    PDV-GL-as-a-model are deliberately not run (see quick_calibrate).
    """
    if H_grid   is None: H_grid   = TRUE_ALPHA_R_GRID
    if eta_grid is None: eta_grid = [0]
    if eta_fixed is None: eta_fixed = eta_grid[0]

    global _ESN_PARAMS
    _ESN_PARAMS = _build_esn(_ARCH)
    assert _ESN_PARAMS["k0"] < 0
    print(f"ESN kappa0 = {_ESN_PARAMS['k0']:.6f}  ✓  (negative leverage guaranteed)\n")

    results = {}

    for eta in eta_grid:          # unused (single dummy value)
        for H in H_grid:          # H == alpha_r, the DGP's true rough exponent
            key = (H, eta)
            print(f"{'='*60}  PDV-GL Table 8, Parameter Set 1  "
                  f"(beta0={TRUE_4F_BETA0}, beta1={TRUE_4F_BETA1}, "
                  f"beta2={TRUE_4F_BETA2}, lam1=[{TRUE_4F_LAM10},{TRUE_4F_LAM11}], "
                  f"theta1={TRUE_4F_THETA1}, lam2=[{TRUE_4F_LAM20},{TRUE_4F_LAM21}], "
                  f"theta2={TRUE_4F_THETA2})")

            # ── DGP (PDV-GL, ground truth) ─────────────────────────
            print("  DGP (PDV-GL, ground truth) ...")
            dgp_sts = dgp_pdvgl(H, n_dgp, T, burn, dt=dt,
                                 seed_base=int(H*1000))
            def _m(k): return float(np.nanmean([s[k] for s in dgp_sts]))
            print(f"  DGP: H_hat={_m('H_hat'):.4f}  vol={_m('mean_vol_ann')*100:.1f}%  "
                  f"Z={_m('zumbach'):.4f}  L={_m('leverage'):.4f}  K={_m('kurtosis'):.3f}")

            # ── Calibrate (ESN + QRH) ───────────────────────────────
            print("  Calibrating ESN + QRH (full 11-term smooth score S) ...")
            esn_cal, qrh_cal, ref = quick_calibrate(dgp_sts, n_cal, T_cal, burn_cal, dt)

            # ── Simulate final paths ──────────────────────────────
            print("  Simulating final paths ...")
            esn_sts = sim_model_paths("ESN A2-981003", esn_cal, n_sim, T, burn, dt)
            qrh_sts = sim_model_paths("QRH", qrh_cal, n_sim, T, burn, dt)

            results[key] = {
                "DGP (PDV-GL)":  dgp_sts,
                "ESN A2-981003": esn_sts,
                "QRH":           qrh_sts,
                "ref":           ref,
            }

            # ── Summary ───────────────────────────────────────────
            print(f"\n  Summary  true_alpha_r={H}:")
            print(f"  {'Model':<22} {'H_hat':>7} {'vol%':>6} {'Z':>8} "
                  f"{'Lev':>8} {'Kurt':>7} {'Score':>7}")
            print("  "+"-"*65)
            for name in ["DGP (PDV-GL)"]+MODEL_NAMES:
                sts=results[key][name]
                sc=float(np.nanmean([score_fn(s,ref) for s in sts]))
                print(f"  {name:<22} {_m('H_hat'):>7.4f} "
                      f"{float(np.nanmean([s['mean_vol_ann'] for s in sts]))*100:>5.1f}% "
                      f"{float(np.nanmean([s['zumbach'] for s in sts])):>8.4f} "
                      f"{float(np.nanmean([s['leverage'] for s in sts])):>8.4f} "
                      f"{float(np.nanmean([s['kurtosis'] for s in sts])):>7.3f} {sc:>7.3f}")
            print()

    # ── Figures ──────────────────────────────────────────────────────────
    print(f"\n{'='*60}\nProducing figures ...\n{'='*60}")

    for eta in eta_grid:
        H_ok = [H for H in H_grid if (H,eta) in results]
        if not H_ok: continue
        sfx = f"eta{str(eta).replace('.','p')}"

        fig = plot_Z1_curves(results, H_ok, eta,
              save_path=f"{save_prefix}_Z1_{sfx}.png")
        plt.close(fig)

        fig = plot_Z3_diagnostics(results, H_ok, eta,
              save_path=f"{save_prefix}_Z3_{sfx}.png")
        plt.close(fig)

        fig = plot_zumbach_figure(results, H_ok, eta,
              save_path=f"{save_prefix}_zumbach_{sfx}.png")
        plt.close(fig)

        fig = plot_leverage_figure(results, H_ok, eta,
              save_path=f"{save_prefix}_leverage_{sfx}.png")
        plt.close(fig)

    fig = plot_Z2_bars(results, H_grid, eta_grid,
          save_path=f"{save_prefix}_Z2.png")
    plt.close(fig)

    fig = plot_Z4_summary(results, H_grid, eta_grid,
          save_path=f"{save_prefix}_Z4.png")
    plt.close(fig)

    # ── LaTeX / PDF ───────────────────────────────────────────────────────
    lp = latex_path or f"{save_prefix}_protocol.tex"
    write_latex(lp)
    if os.system("which pdflatex>/dev/null 2>&1") == 0:
        os.system(f"pdflatex -interaction=nonstopmode {lp}>/dev/null 2>&1")
        os.system(f"pdflatex -interaction=nonstopmode {lp}>/dev/null 2>&1")
        pdf = lp.replace(".tex",".pdf")
        if os.path.exists(pdf): print(f"  Compiled -> {pdf}")

    return results


# ============================================================
# 9.  CLI
# ============================================================

if __name__ == "__main__":
    import argparse
    pa = argparse.ArgumentParser(
        description="Zumbach effect under the PDV-GL model used as ground-truth DGP",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pa.add_argument("--H",   nargs="+", type=float, default=[1],
                    help="PDV-GL true alpha_r (short-memory/rough exponent) value(s)")
    pa.add_argument("--eta", nargs="+", type=float, default=[0],
                    help="unused (kept for CLI compatibility)")
    pa.add_argument("--rho", type=float, default=None,
                    help="unused (PDV-GL DGP leverage comes from its own F(alpha_r) factor)")
    pa.add_argument("--n_dgp", type=int, default=30)
    pa.add_argument("--n_cal", type=int, default=5)
    pa.add_argument("--n_sim", type=int, default=30)
    pa.add_argument("--T",   type=int,   default=2500)
    pa.add_argument("--burn",type=int,   default=500)
    pa.add_argument("--dt",  type=float, default=1.0)
    pa.add_argument("--out", type=str,   default="zumbach_pdvgl")
    pa.add_argument("--fast",action="store_true",
                    help="alpha_r=[0.05,0.15], n=8, T=800")
    args = pa.parse_args()
    if args.fast:
        args.H=[1]; args.eta=[0]
        args.n_dgp=8; args.n_cal=3; args.n_sim=8
        args.T=800; args.burn=200
    run(H_grid=args.H, eta_grid=args.eta, rho=args.rho,
        n_dgp=args.n_dgp, n_sim=args.n_sim, T=args.T, burn=args.burn,
        n_cal=args.n_cal, T_cal=min(args.T,600), burn_cal=min(args.burn,150),
        dt=args.dt, save_prefix=args.out)
