"""Design audit only: a bounded variance correction preserving quoted VIX states.

For the unchanged Gaussian factor generator L and smooth g, h=(d_t+L)g
satisfies E_t int h = E_t g(end)-g(start). Setting g to zero at all quoted
VIX window endpoints protects the model's VIX values at those dates.
This is our proposed M2 extension, not the VDV paper's implementation.
"""
import hashlib,json
from pathlib import Path
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.integrate import quad
from .m2_spx import frozen_stage3
from .full_readout_fit import HERE,dump
from esn_eq.pricing import DatedM2Pricer
OUT=HERE/'runs/m2_protected_variance_design_001'

def run():
    OUT.mkdir(exist_ok=False);architecture,_,premium,_=frozen_stage3();dated=DatedM2Pricer(architecture,premium);kappa=float(dated.layout.rates[6]);times=np.array([27,55,118,153])/365;D=30/365;ends=np.unique(np.r_[times,times+D]);tau=10/365
    # If f_j=tanh((t-b_j)/tau)^2, |f_j'| <= (1-f_j^2)/tau.
    # Products retain this bound: (1-x^2)y+x(1-y^2) <= 1-x^2*y^2,
    # since the difference is (1-x)(1-y)(1-xy)>=0 on[0,1]^2.
    # Thus |a'|<=1/tau independently of the number of endpoints.
    # For phi=tanh(u), |u|sech^2(u)<=1/2 and
    # |tanh(u)|sech^2(u)<=2/(3sqrt3), so |Lphi|<=kappa*(1/2+4/(3sqrt3)).
    K=1/tau+kappa*(.5+4/(3*np.sqrt(3)))
    def envelope(t):
        t=np.asarray(t);z=np.tanh((t[...,None]-ends)/tau);f=z*z;der=2*z*(1-z*z)/tau;a=np.prod(f,axis=-1);ap=sum(der[...,j]*np.prod(np.delete(f,j,axis=-1),axis=-1) for j in range(len(ends)));return a,ap
    def h(t,u):
        a,ap=envelope(t);z=np.tanh(u);d=1-z*z;return (ap*z+a*(-kappa*u*d-2*kappa*z*d))/K
    nodes={n:hermgauss(n) for n in [64,96]}
    def conditional(s,t,u,n):
        z,w=nodes[n];lag=s-t;state=np.exp(-kappa*lag)*u+np.sqrt(-np.expm1(-2*kappa*lag))*np.sqrt(2)*z;return float(np.dot(w/np.sqrt(np.pi),h(s,state)))
    records=[]
    for n in nodes:
        for t in times:
            for u in [-4.,-1.,0.,1.,4.]:
                value,error=quad(lambda s:conditional(s,t,u,n),t,t+D,epsabs=1e-11,epsrel=1e-10,points=ends[(ends>t)&(ends<t+D)].tolist(),limit=100);records.append(dict(GH=n,day=float(t*365),centered_factor=u,conditional_window_integral=value,estimated_quad_error=error))
    residual=max(abs(r['conditional_window_integral']) for r in records);assert residual<1e-9
    t=np.linspace(0,1030/365,1001);u=np.linspace(-8,8,201);boundprobe=float(abs(h(t[:,None],u[None,:])).max());assert 0<boundprobe<=1
    z,w=nodes[96];meanprobe=max(abs(float(np.dot(w/np.sqrt(np.pi),h(tt,np.sqrt(-np.expm1(-2*kappa*tt))*np.sqrt(2)*z)))) for tt in t[::20]);assert meanprobe<1e-14
    off=quad(lambda s:conditional(s,40/365,-1.,96),40/365,70/365,epsabs=1e-10)[0];assert abs(off)>1e-6
    record=dict(status='completed algebraic design audit, no option fit or model adoption',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),baseline='unchanged protected two-speed M2',factor_mode_index=6,half_life_days=float(np.log(2)*365/kappa),kappa=kappa,quoted_VIX_days=(times*365).tolist(),window_endpoints_days=(ends*365).tolist(),envelope_width_days=10,analytic_global_absolute_h_bound=1.,normalization_K=K,observed_grid_max_abs_h=boundprobe,max_conditional_window_integral_error=residual,max_unconditional_mean_h_error=meanprobe,off_quote_window_integral_at_day40_state_minus1=float(off),records=records,positivity='v_new=v_base+epsilon*h, |h|<=1 and v_base>=f>0. If |epsilon|<=f/2 then v_new>=f/2. The original exact hard floor f is not automatically retained.',forward_variance='centered Gaussian mode plus odd tanh and odd OU-generator term gives E[h_t]=0, preserving E[v_t]=xi(t) in the ideal normalized baseline.',protection_scope='model-derived VIX value at four quoted expiries for the same full Gaussian filtration, hence all strikes and futures at those expiries; not unquoted expiries or full VIX process law',calendar_constraint='endpoint zeros are calibration-grid constraints, not a special rule for one maturity; broader expiry coverage needs new constraints or a more general marginal-preserving calibration',evidence_scope='Dynkin argument plus finite Gaussian/quadrature checks; no new Lean proof, compiled pricing correspondence, population-forward certificate or SPX improvement claimed',adopted=False)
    dump(OUT/'design_audit.json',record);(OUT/'source.py').write_text(Path(__file__).read_text());print('PROTECTED VARIANCE DESIGN PASS',residual,boundprobe,meanprobe,'offquoted',off,flush=True)
if __name__=='__main__':run()
