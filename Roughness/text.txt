"""
rough_bergomi_zumbach.py
========================

DGP: Rough Bergomi model (Bayer, Friz & Gatheral 2016).

    dS_t / S_t = sqrt(V_t) dW_t^(1)

    V_t = xi_0 * exp( eta * W_t^H  -  eta^2 / 2 * t^{2H} )

where W_t^H is the Riemann-Liouville (RL) fractional process:

    W_t^H = sqrt(2H) * int_0^t (t-s)^{H-1/2} dW_s^(2)

and  d<W^(1), W^(2)>_t = rho dt  (leverage correlation).

The centering term  -eta^2/2 * t^{2H}  ensures E[V_t] = xi_0.

Parameters explored:
    H   in {0.05, 0.10, 0.20}   rough regime
    eta in {0.5, 1.0, 2.0}      vol-of-vol
    rho = -0.70  (fixed)         spot-vol leverage

Analytical Zumbach property of rBergomi:
    Z(L) = Corr(past_R^2, future_V) - Corr(past_V, future_R^2) > 0
    because rho < 0 and rough memory => past large shocks predict
    future vol (Zumbach effect). The statistic is most pronounced for
    small H (rough) and large eta.

THREE MODELS (same as roughness_evaluation.py — UNCHANGED):
    (A) ESN A2-981003   — notebook-exact, Gaussian innovations
    (B) QRH             — Euler-Volterra, ring buffer 252 steps
    (C) PDV-GL          — two-factor power-law kernel

CALIBRATION (UNCHANGED):
    Full 11-term data-adaptive score S (Eq. 17):
    S = 2.2f(H) + 1.4f(vol) + 0.6g(q995) + 0.2g(maxV)
      + 0.8f(V_ACF) + 1.0f(G_T) + 0.8f(F_T) + 1.0f(Z)
      + 1.1f(L) + 0.7f(K) + 0.5f(A) - 5*Stress
    All centres and tolerances are data-adaptive (DGP cross-path mean/std).

ZUMBACH FIGURES:
  Figure Z1  — Z(L) curves: E[past_R^2 * fut_V] - E[past_V * fut_R^2]
               vs lag L, one column per H, rows = DGP + 3 models.
               Includes ±1 std band across paths and a smoothed mean line.

  Figure Z2  — Zumbach scalar bar chart:
               Z = mean(Z(5), Z(10), Z(20)) per model, one panel per H,
               grouped bars across eta values.

  Figure Z3  — Return / vol diagnostics (histogram, Hill, leverage ACF,
               volatility ACF) — confirming DGP realism.

  Figure Z4  — Summary: Zumbach scalar vs H (fixed eta), vs eta (fixed H),
               score S vs H.

OUTPUT FILES (prefix roughness_bergomi):
  rough_bergomi_zumbach_Z1_{eta}.png
  rough_bergomi_zumbach_Z2.png
  rough_bergomi_zumbach_Z3_{eta}.png
  rough_bergomi_zumbach_Z4.png
  rough_bergomi_zumbach_protocol.tex / .pdf
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

H_GRID   = [0.05, 0.10, 0.20]
ETA_GRID = [0.5, 1.0, 2.0]
RHO_DGP  = -0.70          # fixed leverage

COLORS = {
    "DGP (rBergomi)": "#888780",
    "ESN A2-981003":  "#1a4d80",
    "QRH":            "#D85A30",
    "PDV-GL":         "#1D9E75",
}
CORAL = "#c0392b"; GRAY = "#888780"
MODEL_NAMES = ["ESN A2-981003", "QRH", "PDV-GL"]

_ARCH = dict(
    matrix_seed=202695547565, n_r=64, n_z=12, H_target=0.08,
    rough_scale=0.40, z_strength=0.34, z_readout=0.05,
    even_strength=1.50, linear_strength=0.25, gamma_norm=1.00,
    local_z_strength=0.03, zz_scale=0.08,
    sign_prob_neg=0.22, rough_orientation=-1.0,
)

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
        logret=x, sigma_daily=sigma_daily,
        vol_acf_vec=vol_acf, lev_vec=lev,
    )

def make_score_ref(dgp_stats):
    """Data-adaptive score reference — identical to roughness_evaluation.py."""
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

def score_fn(st, ref):
    """Full 11-term score S — identical to roughness_evaluation.py."""
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
# 2.  DGP — Rough Bergomi
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


def dgp_rough_bergomi(H, eta, rho, n_paths, T, burn,
                      xi0=None, dt=1.0, seed_base=0):
    """
    Simulate the rough Bergomi DGP.

    Model (Bayer, Friz & Gatheral 2016):
        dS_t / S_t = sqrt(V_t) dW_t^(1)
        V_t = xi_0 * exp( eta * W_t^H  -  eta^2/2 * t^{2H} )
        W_t^H = sqrt(2H) * int_0^t (t-s)^{H-1/2} dW_s^(2)
        d<W^(1), W^(2)>_t = rho dt

    Discretisation:
        W^H is drawn via Cholesky of the fBM covariance (exact same exponent H).
        Correlated Brownians:
            dW^(1) = eps1 * sqrt(dt)
            dW^(2) = (rho*eps1 + sqrt(1-rho^2)*eps2) * sqrt(dt)
        Log-variance: log(V_t) = log(xi0) + eta*W^H_t - eta^2/2 * t^{2H}
        Log-return:   dX_t = -V_t/2 * dt + sqrt(V_t) * sqrt(dt) * eps1

    Note on sign: dS/S = sqrt(V) dW^(1), so for the log-return:
        dX = d log S = -V/2 dt + sqrt(V) dW^(1)   (Ito)
    The -V/2 Ito correction IS required here because we track log S, not S.

    Parameters
    ----------
    H        : Hurst exponent of W^H
    eta      : vol-of-vol (> 0)
    rho      : spot-vol correlation (< 0 for leverage)
    n_paths  : Monte Carlo paths
    T        : total trading days (including burn)
    burn     : burn-in days discarded
    xi0      : initial variance (default TV_DAY^2)
    dt       : time step (1 = daily)
    seed_base: base RNG seed

    Returns
    -------
    st_all : list of stat dicts per path (post-burn), each containing all 11
             score statistics PLUS logret, sigma_daily, vol_acf_vec, lev_vec.
    """
    if xi0 is None: xi0=TV_DAY**2

    L=_get_rl_chol(T,H,dt)
    sqrt1mrho2=math.sqrt(max(1.-rho**2,0.))

    rng_base=np.random.default_rng(seed_base); st_all=[]

    for p in range(n_paths):
        rng=np.random.default_rng(int(rng_base.integers(1<<31)))

        eps1=rng.standard_normal(T)   # spot Brownian increments
        eps2=rng.standard_normal(T)   # independent vol Brownian increments

        # Correlated vol driver: dW^(2) = rho*eps1 + sqrt(1-rho^2)*eps2
        dW2=rho*eps1+sqrt1mrho2*eps2

        # RL process W^H_t via fBM Cholesky on the vol driver
        WH=L@dW2           # (T,)  values of W^H at t=dt,2dt,...,T*dt

        # Log-variance:  log V_t = log(xi0) + eta*W^H_t - eta^2/2 * t^{2H}
        t_vec=np.arange(1,T+1,dtype=float)*dt
        log_V=math.log(xi0)+eta*WH-0.5*eta**2*t_vec**(2*H)
        V=np.exp(log_V)                       # (T,) daily variance

        # Log-returns: dX = -V/2 dt + sqrt(V) sqrt(dt) eps1  (Ito on log S)
        sdt=math.sqrt(dt)
        X= -0.5*V*dt + np.sqrt(V)*sdt*eps1   # (T,)

        xb=X[burn:]; vb=V[burn:]
        st=compute_statistics(xb,vb)
        st_all.append(st)

    return st_all

# ============================================================
# 3.  Models — UNCHANGED from roughness_evaluation.py
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
    """Simulate one model and compute all 11 score statistics."""
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
    Calibrate ESN, QRH, PDV-GL to the rough Bergomi DGP.
    Uses the FULL 11-term data-adaptive score S (Eq. 17).
    All calibration logic is UNCHANGED from roughness_evaluation.py.
    """
    global _ESN_PARAMS
    ref=make_score_ref(dgp_sts); tH=ref["H_hat_c"]
    print(f"    Score ref: H={tH:.4f}  vol={ref['mean_vol_ann_c']*100:.1f}%  "
          f"Z={ref['zumbach_c']:.4f}  L={ref['leverage_c']:.4f}  K={ref['kurtosis_c']:.3f}")

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

# ============================================================
# 6.  Figures
# ============================================================

ALL_NAMES=["DGP (rBergomi)"]+MODEL_NAMES
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
        f"\nRough Bergomi DGP  ($\\eta={eta_fixed}$,  $\\rho={RHO_DGP}$)"
        r"  |  $Z(L)>0$ = Zumbach effect present",
        fontsize=11,fontweight="bold")

    for ci,H in enumerate(Hcols):
        data=results[(H,eta_fixed)]
        axes[0,ci].set_title(f"$H={H}$",fontsize=11,fontweight="bold")
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
            ax.set_title("" if ri>0 else f"$H={H}$",fontsize=10 if ri==0 else 9)
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
        "— Rough Bergomi DGP",fontsize=11,fontweight="bold")

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
    ax.set_title(f"$Z$ vs $H$  (fixed $\\eta={eta_fixed}$)",fontsize=10)
    ax.legend(fontsize=9); ax.grid(True,alpha=0.22)
    for H in H_grid:
        if (H,eta_fixed) in results:
            dgpZ=float(np.nanmean([s["zumbach"] for s in results[(H,eta_fixed)]["DGP (rBergomi)"]]))
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
    ax.set_xlabel("Vol-of-vol $\\eta$",fontsize=10)
    ax.set_ylabel("Zumbach scalar $Z$",fontsize=10)
    ax.set_title(f"$Z$ vs $\\eta$  (fixed $H={H_fixed}$)",fontsize=10)
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
        f"DGP + Model Diagnostics — Rough Bergomi  ($\\eta={eta_fixed}$,  $\\rho={RHO_DGP}$)\n"
        "Histogram | Hill | Leverage ACF | Volatility ACF",
        fontsize=11,fontweight="bold")

    for ci,H in enumerate(Hcols):
        data=results[(H,eta_fixed)]
        axes[0,ci].set_title(f"$H={H}$",fontsize=11,fontweight="bold")

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
        "Summary — Rough Bergomi DGP"
        r"  |  $V_t=\xi_0\exp(\eta W_t^H-\eta^2t^{2H}/2)$",
        fontsize=11,fontweight="bold")

    for ci,eta in enumerate(cols_eta):
        Hok=[H for H in H_grid if (H,eta) in results]
        if not Hok: continue

        axes[0,ci].set_title(f"$\\eta={eta}$",fontsize=11,fontweight="bold")

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
{\large Zumbach Effect under the Rough Bergomi DGP}\\[0.2cm]
{\normalsize Othmane Zarhali --- Paris Dauphine / CNRS}
\end{center}
\vspace{0.15cm}\noindent\rule{\textwidth}{0.8pt}

\begin{abstract}
We evaluate the Zumbach time-reversal asymmetry effect produced by three
stochastic volatility models under the rough Bergomi data-generating process (DGP).
The DGP is parameterised by Hurst exponent $H\in\{0.05,0.10,0.20\}$ and
vol-of-vol $\eta\in\{0.5,1.0,2.0\}$ with fixed leverage $\rho=-0.70$.
The three models --- ESN A2-981003, QRH, and PDV-GL --- are calibrated using
the identical full 11-term data-adaptive score $S$ (Eq.~\eqref{eq:score}).
The Zumbach effect is assessed via the per-lag curve $Z(L)$ and the scalar
$Z=\frac{1}{3}[Z(5)+Z(10)+Z(20)]$, supplemented by leverage ACF, vol ACF,
and return tail diagnostics.
\end{abstract}

\tableofcontents\vspace{0.4cm}

\section{DGP: Rough Bergomi Model}

The data-generating process is the rough Bergomi model
\cite{bayer2016,gatheral2018}:
\begin{equation}
  \frac{dS_t}{S_t} = \sqrt{V_t}\,dW_t^{(1)},
  \label{eq:spot}
\end{equation}
\begin{equation}
  V_t = \xi_0\,\exp\!\Bigl(\eta\,\widehat W_t^H
        - \frac{\eta^2}{2}\,t^{2H}\Bigr),
  \label{eq:var}
\end{equation}
where $\widehat W_t^H$ is the Riemann-Liouville (RL) process:
\begin{equation}
  \widehat W_t^H = \sqrt{2H}\int_0^t(t-s)^{H-1/2}\,dW_s^{(2)},
  \label{eq:rl}
\end{equation}
$(W^{(1)},W^{(2)})$ is a 2-D Brownian motion with
$d\langle W^{(1)},W^{(2)}\rangle_t=\rho\,dt$, $\rho<0$.

\paragraph{Centering.}
The term $-\eta^2 t^{2H}/2$ ensures $\mathbb{E}[V_t]=\xi_0$ for all $t$,
since $\mathbb{E}[\exp(\eta\widehat W_t^H)]=\exp(\eta^2\mathrm{Var}(\widehat W_t^H)/2)
=\exp(\eta^2 t^{2H}/2)$.

\paragraph{Log-return (Ito's lemma).}
Since $dS_t/S_t=\sqrt{V_t}\,dW_t^{(1)}$, the Ito formula gives:
\begin{equation}
  dX_t = d\log S_t = -\tfrac{1}{2}V_t\,dt + \sqrt{V_t}\,dW_t^{(1)}.
  \label{eq:logret}
\end{equation}
The $-V_t/2$ correction is mandatory; it is \emph{not} a modelling choice.

\paragraph{Zumbach effect mechanism.}
The RL rough kernel gives $V_t$ a long memory of past shocks.
Because $\rho<0$, a large negative spot return $dW^{(1)}<0$ drives $V$ upward
(through the correlated $dW^{(2)}=\rho\,dW^{(1)}+\cdots$).
This asymmetry makes past squared returns predictive of future variance but
not vice versa, producing $Z(L)>0$.
The roughness $H<\tfrac{1}{2}$ amplifies this at short lags.

\paragraph{Simulation.}
The correlated Brownians are generated as:
\begin{equation}
  dW^{(1)} = \varepsilon_1\sqrt{dt},\quad
  dW^{(2)} = \bigl(\rho\varepsilon_1+\sqrt{1-\rho^2}\,\varepsilon_2\bigr)\sqrt{dt},\quad
  \varepsilon_1,\varepsilon_2\overset{\mathrm{iid}}{\sim}\mathcal{N}(0,1).
\end{equation}
$\widehat W^H$ is simulated via Cholesky of the fBM covariance matrix
$\Sigma_{ij}=\tfrac{1}{2}(t_i^{2H}+t_j^{2H}-|t_i-t_j|^{2H})$,
which has the same H\"older exponent $H$ as the RL process.

\paragraph{Parameter grid.}
$H\in\{0.05,0.10,0.20\}$,\quad
$\eta\in\{0.5,1.0,2.0\}$,\quad
$\rho=-0.70$ (fixed),\quad
$\xi_0=(0.20/\sqrt{252})^2$.

\section{Models (Unchanged)}

All three models use Gaussian Brownian innovations.

\subsection{ESN A2-981003}
Exact notebook simulate\_path loop ($N_r=64$, $N_z=12$,
$\kappa_0=-0.472416<0$): single $\varepsilon_t\sim\mathcal{N}(0,1)$ drives both
the Ito log-return and the diagonal OU rough bank.
The negative instantaneous leverage $\kappa_0$ gives the ESN an endogenous
Zumbach-like asymmetry.

\subsection{Quadratic Rough Heston (QRH)}
Bourgey \& Gatheral \cite{bourgey2026}: $V_t=Y_t^2+c$, gamma kernel
$\kappa(\tau)=\nu/\Gamma(\alpha)\tau^{\alpha-1}e^{-\lambda\tau}$,
$\alpha=H+\tfrac{1}{2}$, Euler-Volterra with 252-step ring buffer.

\subsection{PDV Guyon-Lekeufack (GL)}
Guyon \& Lekeufack \cite{guyon2023}: two-factor power-law kernel
$V_t=(1-\beta)\bar V+\beta[\alpha_1 F_t(\alpha_r)+(1-\alpha_1)F_t(\alpha_v)]$,
online from full return history.

\section{Calibration (Unchanged --- Full 11-Term Score)}

The data-adaptive score (Eq.~17 of the protocol):
\begin{equation}
\begin{aligned}
S &= 2.2\,f(\hat H)+1.4\,f(\bar\sigma)+0.6\,g(q_{995})+0.2\,g(\sigma_{\max})\\
  &+0.8\,f(\bar V_{\mathrm{ACF}})+1.0\,f(G_T)+0.8\,f(F_T)+1.0\,f(Z)\\
  &+1.1\,f(L)+0.7\,f(K)+0.5\,f(A)-5\cdot\mathrm{Stress},
\end{aligned}
\label{eq:score}
\end{equation}
where $f(x,c,s)=\max(0,1-|x-c|/s)$ and $g(x,c,s)=\max(0,1-(x-c)_+/s)$.
All centres $c_k$ and tolerances $s_k$ are DGP cross-path mean and
$\max(\mathrm{std},0.30|c_k|,\mathrm{floor}_k)$.
The Zumbach term $f(Z,c_Z,s_Z)$ (weight 1.0) penalises any model that
fails to reproduce the DGP's Zumbach magnitude and sign simultaneously.

Calibrated parameters:
\begin{center}\renewcommand{\arraystretch}{1.4}
\begin{tabular}{lll}\toprule
Model & Parameters \\\midrule
ESN A2-981003 & $\Delta b_0$, vol scale $s$ \\
QRH & $H,\nu_{\rm vol},\lambda,c_{\rm frac}$ \\
PDV-GL & $\beta,\alpha_1,\alpha_r,\alpha_v$ \\
\bottomrule\end{tabular}\end{center}

\section{Zumbach Diagnostic}

Two complementary statistics are computed.

\subsection{Windowed integrated-variance scalar (calibration score)}

\begin{equation}
  Z(L) = \mathrm{Corr}\!\bigl(\mathrm{past}_L R^2,\,\mathrm{fut}_L V\bigr)
        - \mathrm{Corr}\!\bigl(\mathrm{past}_L V,\,\mathrm{fut}_L R^2\bigr)>0,
  \label{eq:zumbach_scalar}
\end{equation}
where $\mathrm{past}_L R^2=(X_t{-}X_{t-L})^2$ and
$\mathrm{fut}_L V=(V_{t+L}{-}V_t)/L$.
Scalar $Z=\tfrac{1}{3}[Z(5){+}Z(10){+}Z(20)]$ enters the calibration score
(weight 1.0).

\subsection{Cross-moment statistic (main figure)}

\begin{equation}
  Z^{\rm cm}(\tau)
  = \mathbb{E}\bigl[r^2_{t+\tau}\,r_t\bigr]
  - \mathbb{E}\bigl[r_{t+\tau}\,r_t^2\bigr] < 0,
  \label{eq:zumbach_cm}
\end{equation}
with $\sigma^3$-normalised decomposition:
\begin{equation}
  Z_+(\tau) = \mathbb{E}[r^2_{t+\tau}\,r_t]/\sigma^3, \qquad
  Z_-(\tau) = \mathbb{E}[r_{t+\tau}\,r_t^2]/\sigma^3.
\end{equation}
$Z^{\rm cm}(\tau)<0$ signals the Zumbach effect.
With $\rho<0$, large negative $r_t$ predicts high $r^2_{t+\tau}$, making
$Z_+<0$ and $Z^{\rm cm}<0$ at all lags.

\section{Figures}

\textbf{Main figure (Zumbach cross-moment).}
Four rows, one column per $H$ (fixed $\eta$):
\begin{enumerate}[itemsep=2pt]
  \item $Z^{\rm cm}(\tau)$ vs lag: Gaussian-smoothed mean (line,
    $\sigma_{\rm sm}=2$d) + dots + $\pm1$std band. Negative = Zumbach.
  \item Grouped bars of $Z^{\rm cm}(\tau)$ at lags 1d, 5d, 10d, 20d.
  \item Decomposition: solid $Z_+(\tau)/\sigma^3$, dashed $Z_-(\tau)/\sigma^3$.
  \item Horizontal bar of $Z^{\rm cm}(1\mathrm{d})$, ranked most negative at
    bottom; values annotated in scientific notation.
\end{enumerate}

\textbf{Figure Z1.} Windowed $Z(L)$ curves (row per source, column per $H$).

\textbf{Figure Z2.} Scalar $Z$ vs $H$ and vs $\eta$ for all sources.

\textbf{Figure Z3.} Return diagnostics: histogram | Hill | leverage ACF | vol ACF.

\textbf{Figure Z4.} Summary: $Z$, leverage, $\hat H$, score $S$ vs $H$.

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
$H$ grid & $\{0.05,0.10,0.20\}$ \\
$\eta$ grid & $\{0.5,1.0,2.0\}$ \\
$\rho$ & $-0.70$ (fixed) \\
Smoothing $\sigma$ & 2 days \\
\bottomrule\end{tabular}\end{center}

\begin{thebibliography}{9}
\bibitem{bayer2016} C.\ Bayer, P.\ Friz, J.\ Gatheral.
\textit{Pricing under rough volatility}. QF 16(6):887--904, 2016.
\bibitem{gatheral2018} J.\ Gatheral, T.\ Jaisson, M.\ Rosenbaum.
\textit{Volatility is rough}. QF 18(6):933--949, 2018.
\bibitem{bourgey2026} F.\ Bourgey, J.\ Gatheral.
\textit{Quadratic Rough Heston}. SSRN:5239929, 2026.
\bibitem{guyon2023} J.\ Guyon, J.\ Lekeufack.
\textit{Volatility is (mostly) path-dependent}. QF 23(9), 2023.
\bibitem{zumbach2009} G.\ Zumbach.
\textit{Time reversal invariance in finance}. QF 9(5):505--515, 2009.
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
    Exact replica of the uploaded figure layout.

    One column per H value (fixed eta).  Four rows:

    Row 0 — Z(tau) curves
        Per-lag cross-moment Z(tau) = E[r^2_{t+tau} r_t] - E[r_{t+tau} r^2_t].
        Gaussian-smoothed mean line + unsmoothed dots + ±1std shaded band.
        Reference at Z=0. Negative = Zumbach effect.
        Subtitle: "Z(tau) curves / (most negative at short lag)"

    Row 1 — Z at selected lags  (grouped bar chart)
        Bars at lags 1d, 5d, 10d, 20d (smoothed values).
        One bar per source per lag group.

    Row 2 — Decomposition
        Solid lines:  E[r^2_{t+tau} r_t] / sigma^3  (Zpos)
        Dashed lines: E[r_{t+tau} r^2_t] / sigma^3  (Zneg)
        Subtitle: "Decomp.: solid E[r^2 r']/sigma^3, dashed E[r r'^2]/sigma^3"

    Row 3 — Z(1d) ranked horizontal bar
        Single bar per source showing smoothed Z at lag=1.
        Sorted so most negative is at bottom (strongest Zumbach).
        Annotation with numeric value.
    """
    from scipy.ndimage import gaussian_filter1d

    Hcols = [H for H in H_grid if (H, eta_fixed) in results]
    n_c   = len(Hcols)
    lags  = np.arange(1, max_lag + 1)

    fig, axes = plt.subplots(4, n_c, figsize=(6*n_c, 18))
    if n_c == 1: axes = axes[:, None]

    fig.suptitle(
        f"Zumbach Time-Reversal Asymmetry — Rough Bergomi DGP"
        f"  ($\\eta={eta_fixed}$,  $\\rho={RHO_DGP}$)\n"
        r"$Z(\tau) = E[r^2_{t+\tau}\,r_t] - E[r_{t+\tau}\,r^2_t]$"
        r"  —  $\mathbf{negative}$ for all $\tau$ = Zumbach effect",
        fontsize=12, fontweight="bold")

    for ci, H in enumerate(Hcols):
        data = results[(H, eta_fixed)]

        # Pre-compute curves for all sources
        curves = {}
        for name in ALL_NAMES:
            Zm,Zs,Zpm,Zps,Znm,Zns = zumbach_crossmoment_curve(
                data[name], max_lag)
            # Gaussian smooth
            Zm_sm  = gaussian_filter1d(np.nan_to_num(Zm),  smooth_sigma)
            Zpm_sm = gaussian_filter1d(np.nan_to_num(Zpm), smooth_sigma)
            Znm_sm = gaussian_filter1d(np.nan_to_num(Znm), smooth_sigma)
            # sigma^3 normalisation
            sig3 = float(np.nanmean([
                np.std(s["logret"])**3 for s in data[name]])) + 1e-20
            curves[name] = dict(Zm=Zm, Zs=Zs, Zm_sm=Zm_sm,
                                Zpm_sm=Zpm_sm/sig3, Znm_sm=Znm_sm/sig3)

        col_title = f"$H={H}$,  $\\eta={eta_fixed}$"

        # ── Row 0: Z(tau) curves ──────────────────────────────────────
        ax = axes[0, ci]
        for name in ALL_NAMES:
            c = curves[name]; col = _col(name)
            lw = 2.4 if "DGP" in name else 1.6
            ax.fill_between(lags, c["Zm"]-c["Zs"], c["Zm"]+c["Zs"],
                            alpha=0.13, color=col)
            ax.scatter(lags, c["Zm"], color=col, s=7, alpha=0.35, zorder=3)
            ax.plot(lags, c["Zm_sm"], color=col, lw=lw, label=name, zorder=4)
        ax.axhline(0, color=GRAY, lw=0.9, ls="--", alpha=0.6)
        ax.set_xlabel("Lag (days)", fontsize=9)
        if ci == 0: ax.set_ylabel("$Z(\\tau)$", fontsize=10)
        ax.set_title(f"{col_title}\n$Z(\\tau)$ curves\n(most negative at short lag)",
                     fontsize=9)
        ax.legend(fontsize=7, loc="lower right")
        ax.grid(True, alpha=0.18)

        # ── Row 1: Z at selected lags (grouped bar) ───────────────────
        ax = axes[1, ci]
        n_src = len(ALL_NAMES); bw = 0.70 / n_src
        x_pos = np.arange(len(sel_lags))
        for si, name in enumerate(ALL_NAMES):
            c = curves[name]; col = _col(name)
            vals = [float(c["Zm_sm"][L-1]) for L in sel_lags]
            offset = (si - n_src/2 + 0.5) * bw
            ax.bar(x_pos + offset, vals, bw*0.92,
                   color=col, alpha=0.85, label=name)
        ax.axhline(0, color=GRAY, lw=0.9, ls="--")
        ax.set_xticks(x_pos)
        ax.set_xticklabels([f"{L}d" for L in sel_lags])
        if ci == 0: ax.set_ylabel("$Z(\\tau)$ (smoothed)", fontsize=10)
        ax.set_title("$Z$ at selected lags", fontsize=9)
        ax.legend(fontsize=6); ax.grid(True, axis="y", alpha=0.18)

        # ── Row 2: Decomposition Zpos/sig3 solid, Zneg/sig3 dashed ───
        ax = axes[2, ci]
        for name in ALL_NAMES:
            c = curves[name]; col = _col(name)
            lw = 2.4 if "DGP" in name else 1.5
            ax.plot(lags, c["Zpm_sm"], color=col, lw=lw, ls="-",
                    label=f"{name} $E[r^2r']$")
            ax.plot(lags, c["Znm_sm"], color=col, lw=lw*0.65, ls="--")
        ax.axhline(0, color=GRAY, lw=0.6, ls="--", alpha=0.5)
        ax.set_xlabel("Lag (days)", fontsize=9)
        if ci == 0: ax.set_ylabel("$/\\sigma^3$", fontsize=10)
        ax.set_title("Decomp.: solid $E[r^2r']/\\sigma^3$,"
                     " dashed $E[rr'^2]/\\sigma^3$", fontsize=8)
        ax.legend(fontsize=6); ax.grid(True, alpha=0.18)

        # ── Row 3: Z(1d) ranked horizontal bar ───────────────────────
        ax = axes[3, ci]
        z1d = {name: float(curves[name]["Zm_sm"][0]) for name in ALL_NAMES}
        # sort: most negative first (strongest Zumbach at bottom)
        sorted_items = sorted(z1d.items(), key=lambda x: x[1])
        names_s = [it[0] for it in sorted_items]
        vals_s  = [it[1] for it in sorted_items]
        cols_s  = [_col(n) for n in names_s]
        y_ = np.arange(len(names_s))
        ax.barh(y_, vals_s, color=cols_s, alpha=0.85)
        ax.axvline(0, color=GRAY, lw=0.8, ls="--")
        ax.set_yticks(y_); ax.set_yticklabels(names_s, fontsize=8)
        if ci == 0: ax.set_xlabel("$Z$ at lag = 1 day (smoothed)", fontsize=9)
        ax.set_title("$Z(1d)$ ranked\nmore negative = stronger effect",
                     fontsize=8)
        ax.grid(True, axis="x", alpha=0.18)
        for yi, (name, v) in enumerate(zip(names_s, vals_s)):
            ha = "right" if v < 0 else "left"
            pad = -0.001 if v < 0 else 0.001
            ax.text(v + pad, yi, f"{v:.2e}",
                    ha=ha, va="center", fontsize=7)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
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

    if H_grid   is None: H_grid   = H_GRID
    if eta_grid is None: eta_grid = ETA_GRID
    if eta_fixed is None: eta_fixed = eta_grid[len(eta_grid)//2]

    global _ESN_PARAMS
    _ESN_PARAMS = _build_esn(_ARCH)
    assert _ESN_PARAMS["k0"] < 0
    print(f"ESN kappa0 = {_ESN_PARAMS['k0']:.6f}  ✓  (negative leverage guaranteed)\n")

    results = {}

    for eta in eta_grid:
        for H in H_grid:
            key = (H, eta)
            print(f"{'='*60}  H={H}  eta={eta}  rho={rho}")

            # ── DGP ──────────────────────────────────────────────
            print("  DGP (rough Bergomi) ...")
            dgp_sts = dgp_rough_bergomi(
                H, eta, rho, n_dgp, T, burn,
                dt=dt, seed_base=int(H*1000)+int(eta*100))
            def _m(k): return float(np.nanmean([s[k] for s in dgp_sts]))
            print(f"  DGP: H_hat={_m('H_hat'):.4f}  vol={_m('mean_vol_ann')*100:.1f}%  "
                  f"Z={_m('zumbach'):.4f}  L={_m('leverage'):.4f}  K={_m('kurtosis'):.3f}")

            # ── Calibrate ─────────────────────────────────────────
            print("  Calibrating (full 11-term score S) ...")
            esn_cal, qrh_cal, gl_cal, ref = quick_calibrate(
                dgp_sts, n_cal, T_cal, burn_cal, dt)

            # ── Simulate final paths ──────────────────────────────
            print("  Simulating final paths ...")
            esn_sts = sim_model_paths("ESN A2-981003", esn_cal, n_sim, T, burn, dt)
            qrh_sts = sim_model_paths("QRH",           qrh_cal, n_sim, T, burn, dt)
            gl_sts  = sim_model_paths("PDV-GL",         gl_cal,  n_sim, T, burn, dt)

            results[key] = {
                "DGP (rBergomi)": dgp_sts,
                "ESN A2-981003":  esn_sts,
                "QRH":            qrh_sts,
                "PDV-GL":         gl_sts,
                "ref":            ref,
            }

            # ── Summary ───────────────────────────────────────────
            print(f"\n  Summary  H={H}  eta={eta}:")
            print(f"  {'Model':<22} {'H_hat':>7} {'vol%':>6} {'Z':>8} "
                  f"{'Lev':>8} {'Kurt':>7} {'Score':>7}")
            print("  "+"-"*65)
            for name in ["DGP (rBergomi)"]+MODEL_NAMES:
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
        description="Zumbach effect under rough Bergomi DGP",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pa.add_argument("--H",   nargs="+", type=float, default=[0.05,0.10,0.20])
    pa.add_argument("--eta", nargs="+", type=float, default=[0.5,1.0,2.0])
    pa.add_argument("--rho", type=float, default=-0.70)
    pa.add_argument("--n_dgp", type=int, default=30)
    pa.add_argument("--n_cal", type=int, default=5)
    pa.add_argument("--n_sim", type=int, default=30)
    pa.add_argument("--T",   type=int,   default=2500)
    pa.add_argument("--burn",type=int,   default=500)
    pa.add_argument("--dt",  type=float, default=1.0)
    pa.add_argument("--out", type=str,   default="rough_bergomi_zumbach")
    pa.add_argument("--fast",action="store_true",
                    help="H=[0.05,0.10], eta=[1.0], n=8, T=800")
    args = pa.parse_args()
    if args.fast:
        args.H=[0.05,0.10]; args.eta=[1.0]
        args.n_dgp=8; args.n_cal=3; args.n_sim=8
        args.T=800; args.burn=200
    run(H_grid=args.H, eta_grid=args.eta, rho=args.rho,
        n_dgp=args.n_dgp, n_sim=args.n_sim, T=args.T, burn=args.burn,
        n_cal=args.n_cal, T_cal=min(args.T,600), burn_cal=min(args.burn,150),
        dt=args.dt, save_prefix=args.out)
