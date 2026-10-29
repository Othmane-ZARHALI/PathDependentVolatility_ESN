"""A fixed-marginal mechanism oracle; NOT an M2 calibration or VIX process.

Match E[X]=1 and E[X^2] between a lognormal X and sqrt(two Gaussian
exponentials). Direct quadrature is compared with independent defensive-mixture
Sobol scrambles. Shared scrambles couple the two candidates. No price clipping.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr, logsumexp
from scipy.stats import norm, qmc

HERE = Path(__file__).resolve().parent
A, B, W = .18, 1.8, .1
LIMIT = 12.
T = 153/365
SEEDS = list(range(2026102600, 2026102616))

def response(z):
    return np.sqrt((1-W)*np.exp(A*z-A*A/2)+W*np.exp(B*z-B*B/2))

def tail_bound(cut):
    # sqrt(x+y)<=sqrt(x)+sqrt(y); integrate each branch analytically.
    return sum(np.sqrt(w)*np.exp(-a*a/8)*
               (norm.sf(cut-a/2)+norm.cdf(-cut-a/2))
               for a,w in [(A,1-W),(B,W)])

mean, mean_error = quad(lambda z:response(z)*norm.pdf(z),-LIMIT,LIMIT,epsabs=1e-13,epsrel=1e-13)
second = 1/mean**2
single_sd = np.sqrt(np.log(second))

def black(k, vol_time):
    d1 = -np.log(k)/vol_time+vol_time/2
    return ndtr(d1)-k*ndtr(d1-vol_time)

strikes = np.array([1.,2.,4.,8.])
rows=[]
for k in strikes:
    root=brentq(lambda z:response(z)/mean-k,-LIMIT,LIMIT,xtol=1e-13)
    price,err=quad(lambda z:(response(z)/mean-k)*norm.pdf(z),root,LIMIT,epsabs=1e-14,epsrel=1e-12)
    iv=brentq(lambda v:black(k,v)-price,1e-4,5.,xtol=1e-13)/np.sqrt(T)
    rows.append(dict(K_over_F=float(k),mixed_call=price,quad_reported_error=err,
                     omitted_tail_upper_bound=tail_bound(LIMIT)/mean,
                     mixed_annual_IV=iv,single_call=float(black(k,single_sd)),
                     single_annual_IV=float(single_sd/np.sqrt(T))))

# Exact likelihood ratio against the fixed full-support proposal.
mix=np.array([.4,.2,.2,.2]);shifts=np.array([0.,2.,4.,6.]);blocks=[]
for seed in SEEDS:
    u=qmc.Sobol(2,scramble=True,seed=seed).random_base2(13)
    j=np.searchsorted(np.cumsum(mix),u[:,0],side='right')
    z=norm.ppf(u[:,1])+shifts[j]
    weight=np.exp(norm.logpdf(z)-logsumexp(np.log(mix)[None,:]+norm.logpdf(z[:,None]-shifts),axis=1))
    x=response(z)/mean;y=np.exp(single_sd*z-.5*single_sd**2)
    blocks.append(dict(seed=seed,mass=float(weight.mean()),
        mixed_forward=float(np.mean(weight*x)),mixed_second=float(np.mean(weight*x*x)),
        single_forward=float(np.mean(weight*y)),single_second=float(np.mean(weight*y*y)),
        mixed_calls=np.mean(weight[:,None]*np.maximum(x[:,None]-strikes,0),axis=0).tolist(),
        single_calls=np.mean(weight[:,None]*np.maximum(y[:,None]-strikes,0),axis=0).tolist()))
for i,row in enumerate(rows):
    for name in ['mixed','single']:
        v=np.array([b[name+'_calls'][i] for b in blocks]);se=v.std(ddof=1)/np.sqrt(len(v))
        row[name+'_simulation_call']=float(v.mean());row[name+'_simulation_SE']=float(se)
        row[name+'_simulation_minus_oracle_in_SE']=float((v.mean()-row[name+'_call'])/se) if se>0 else None
        row[name+'_simulation_note']='Zero observed dispersion cannot certify the rare-tail price.' if se==0 else 'Scramble SE is a numerical diagnostic.'

out=dict(scope='Mechanism benchmark only: one marginal, no M2 reservoir, floor, 30-day conditional window, market fitting or time-consistent process.',
    equation='R(z)=sqrt(.9 exp(.18 z-.18^2/2)+.1 exp(1.8 z-1.8^2/2)); X=R/E[R]. Single Y=exp(s z-s^2/2), s^2=log(E[X^2]).',
    formal_correspondence='Lean proves the algebraic log-response curvature numerator and a generic moment obstruction. It does not prove numerical prices, derivatives, Gaussian integrals, implied-volatility slopes or the M2 adaptation.',
    mean_response=mean,mean_quad_reported_error=mean_error,mean_omitted_tail_upper_bound=tail_bound(LIMIT),
    common_forward=1.,common_second_moment=second,illustrative_expiry_days=153,
    numerical_note='QUADPACK errors are diagnostics, not rigorous bounds; omitted-tail expression is an analytic bound evaluated in floating point. Exact Gaussian second moment and forward equalities use exact normalization; implementation uses quadrature mean.',
    mixture_weights=mix.tolist(),mixture_shifts=shifts.tolist(),points_per_scramble=8192,
    results=rows,simulation_blocks=blocks,
    development_note='development_probe.py is exploratory and has an unmatched hard-coded comparator; excluded from final results. A preliminary large-sensitivity trial was rejected because its finite-window calculation did not qualify far-tail contributions.',
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['simulation_blocks']},indent=2))
