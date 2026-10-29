"""Stress-variance-scaled Dynkin correction with the original floor retained."""
from .m2_rich_protected_variance import *
from .m2_rich_protected_variance import Kernel as RichKernel

OUT=HERE/'runs/m2_available_variance_001'
ENDSAV=np.array([27,29,55,57,85,91.25,118,119,148,153,183,211,273.75,302,666,1030])/365
LAM=[.125,.5]

def stress_parameters(x):
    hf=np.exp(np.log(365)*x[21]);hs=hf*(365/hf)**x[18]
    return np.log(2)*365/np.array([hs,hf]),5*x[[16,20]]

def stress_bound(a,b,lam):return a*((4*lam+2)/np.e+b*(1+2*lam)/np.sqrt(2*lam*np.e)+2*lam)

def clock(t,k,D):
    t=np.asarray(t);q=np.zeros_like(t);qd=np.zeros_like(t);B=2*k*np.sqrt(2/(np.pi*np.e))+D
    def put(mask,a,ap,A):
        K=B/2+np.sqrt(A*A+B*B/4);q[mask]=a/K;qd[mask]=ap/K
    m=t<ENDSAV[0];A=np.pi/(2*ENDSAV[0]);put(m,np.cos(A*t[m])**2,-A*np.sin(2*A*t[m]),A)
    for left,right in zip(ENDSAV[:-1],ENDSAV[1:]):
        m=(t>=left)&(t<right);A=np.pi/(right-left);theta=A*(t[m]-left);put(m,np.sin(theta)**2,A*np.sin(2*theta),A)
    m=t>=ENDSAV[-1];T=10/365;A=np.pi/(2*T);u=t[m]-ENDSAV[-1];put(m,np.where(u>=T,1.,np.sin(A*u)**2),np.where(u>=T,0.,A*np.sin(2*u*A)),A)
    return q,qd

def scalar(t,u,z,k,a,b,lam):
    q,qd=clock(t,k,stress_bound(a,b,lam));v=-np.expm1(-2*a*np.asarray(t));y=z-b*v;g=np.exp(-lam*y*y);A=np.exp(b*z-.5*b*b*v);phi=erf(u/np.sqrt(2));lp=-2*k*u*np.sqrt(2/np.pi)*np.exp(-u*u/2);D=a*((4*lam*lam+2*lam)*y*y-b*(1+2*lam*v)*y-2*lam)
    return np.stack(np.broadcast_arrays(q*phi*A*g,A*g*(qd*phi+q*lp+q*phi*D),A),axis=-1)

def conditional(s,T,u,z,k,a,b,lam):
    lag=max(0,s-T);mu=np.exp(-k*lag)*u;vu=-np.expm1(-2*k*lag);ep=erf(mu/np.sqrt(2*(1+vu)));elp=-2*k*np.sqrt(2/np.pi)*mu/(1+vu)**1.5*np.exp(-mu*mu/(2*(1+vu)))
    decay=np.exp(-a*lag);vz=-np.expm1(-2*a*lag);qt=-np.expm1(-2*a*s);qT=-np.expm1(-2*a*T);mean=decay*z-b*decay**2*qT;den=1+2*lam*vz;MA=np.exp(b*decay*z-.5*b*b*decay**2*qT);Eg=np.exp(-lam*mean*mean/den)/np.sqrt(den);m=mean/den;v=vz/den
    EAg=MA*Eg;EAD=EAg*a*((4*lam*lam+2*lam)*(v+m*m)-b*(1+2*lam*qt)*m-2*lam);q,qd=clock(s,k,stress_bound(a,b,lam));return (qd*ep+q*elp)*EAg+q*ep*EAD

class Kernel:
    def __init__(self,path):
        self.rich=RichKernel(path);self.k=self.rich.k;self.k.available_stock.argtypes=self.k.rich_stock.argtypes
        P=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');self.k.available_probe.argtypes=[ct.c_int]+[P]*7
    def probe(self,t,u,z,k,a,b,lam):
        t,u,z=np.broadcast_arrays(t,u,z);out=np.zeros((t.size,3));cfg=arr([k,a,b,lam]);self.k.available_probe(t.size,arr(t.ravel()),arr(u.ravel()),arr(z.ravel()),cfg,cfg,arr(ENDSAV),out);return out.reshape(t.shape+(3,))
    def stocks(self,d,pairs,p,x,configs,rates,seed):
        configs=arr(configs);assert np.max(abs(configs[:,2]))<=1
        models=arr(np.tile(x,(len(configs),1)));out=np.zeros((2*pairs,int(d['expiries']),len(configs)+1,2));self.k.available_stock(len(d['dt']),pairs,d['projection'].shape[1],int(d['expiries']),len(configs),d['transition'],d['decay'],d['root'],d['cross'],d['mean'],d['projection'],d['time'],d['dt'],d['seg'],d['coarse'],d['reference'],d['expiry'],d['lognorm'],arr(p),models,configs,arr(np.r_[rates,ENDSAV]),1.,int(seed),out);assert np.isfinite(out).all() and (out>0).all();return out

def audit(k,d,p,x,rates):
    sr,b=stress_parameters(x);rng=np.random.default_rng(20261029201);records=[];err=bound=fd=oracle=0.
    for rate in rates:
        for lam in LAM:
            for a,bb in zip(sr,b):
                t=rng.uniform(.001,1030/365,100);u=rng.uniform(-4,4,100);z=rng.uniform(-4,7,100);got=k.probe(t,u,z,rate,a,bb,lam);ref=scalar(t,u,z,rate,a,bb,lam);err=max(err,float(np.max(abs(got-ref)/(1+abs(ref)))));bound=max(bound,float(np.max(abs(got[:,1])/got[:,2])))
                # Independently finite-difference the entire state-dependent potential.
                keep=np.min(abs(t[:,None]-np.r_[ENDSAV,ENDSAV[-1]+10/365]),axis=1)>1e-4;t=t[keep];u=u[keep];z=z[keep];dt=1e-7;dx=1e-4
                def pot(t,u,z):return scalar(t,u,z,rate,a,bb,lam)[...,0]
                q=pot(t,u,z);L=(pot(t+dt,u,z)-pot(t-dt,u,z))/(2*dt)
                for state,rr,isU in [(u,rate,True),(z,a,False)]:
                    plus=pot(t,u+dx,z) if isU else pot(t,u,z+dx);minus=pot(t,u-dx,z) if isU else pot(t,u,z-dx);L+=rr*((plus-2*q+minus)/dx**2-state*(plus-minus)/(2*dx))
                actual=k.probe(t,u,z,rate,a,bb,lam);fd=max(fd,float(np.max(abs(L-actual[:,1])/(1+actual[:,2]))))
                for day in [27,55,118,153]:
                    T=day/365
                    for U,Z in [(-3.,-2.),(-1.,0.),(1.,2.),(3.,5.)]:
                        initialA=np.exp(bb*Z-.5*bb*bb*(-np.expm1(-2*a*T)));maxlevel=max(x[:16])
                        def integrand(s):
                            seg=min(np.searchsorted(ENDSAV,s,side='right'),15);return x[seg]*conditional(s,T,U,Z,rate,a,bb,lam)/(initialA*maxlevel)
                        val,e=quad(integrand,T,T+30/365,epsabs=1e-10,points=ENDSAV[(ENDSAV>T+1e-13)&(ENDSAV<T+30/365-1e-13)].tolist(),limit=100);records.append(dict(kappa=float(rate),stress_rate=float(a),lam=lam,day=day,U=U,Z=Z,scaled_integral=val,estimated_quad_error=e))
                # Full two-dimensional Gaussian quadrature versus analytic conditional oracle.
                nodes,weights=hermgauss(64);nodes*=np.sqrt(2);weights/=np.sqrt(np.pi);T=55/365;s=T+.017;lag=s-T;U=-1.;Z=2.;uu=np.exp(-rate*lag)*U+np.sqrt(-np.expm1(-2*rate*lag))*nodes[:,None];zz=np.exp(-a*lag)*Z+np.sqrt(-np.expm1(-2*a*lag))*nodes[None,:];actual=float(weights@k.probe(s,uu,zz,rate,a,bb,lam)[...,1]@weights);expected=float(conditional(s,T,U,Z,rate,a,bb,lam));oracle=max(oracle,abs(actual-expected)/(1+abs(expected)))
    residual=max(abs(r['scaled_integral']) for r in records);assert err<1e-10 and fd<5e-6 and bound<=1 and residual<1e-8 and oracle<1e-8,(err,fd,bound,residual,oracle)
    small=dict(d);n=97
    for key in ['time','dt','mean','seg','coarse','reference','expiry','transition','lognorm']:small[key]=d[key][:n].copy()
    small['expiry'][:]=-1;small['expiry'][-1]=0;small['expiries']=1;got=k.stocks(small,4,p,x,arr([[0,0,0,0],[0,0,0,1]]),rates,20261029202);coeff=np.zeros((2,13));coeff[1,6]=1;old=k.rich.stocks(small,4,p,x,coeff,rates,20261029202);replay=float(abs(got-old).max());assert replay<1e-12
    dump(OUT/'implementation_audit.json',dict(compiled_formula_relative_error=err,finite_difference_full_generator_error_scaled=fd,observed_max_abs_h_over_stress=bound,conditional_scaled_integral_error=residual,conditional_oracle_GH64_error=oracle,baseline_and_rich_replay_error=replay,records=records,analytic_bound='|h_j|<=existing stress variance; |eta|<=1 implies v_new>=original floor f',no_new_Lean_proof=True));print('AVAILABLE AUDIT PASS',err,fd,bound,residual,oracle,replay,flush=True)

def stage(k,d,p,x,rates,configs,name,pairs,tag):
    folder=OUT/name;folder.mkdir();pieces=[];seeds=[]
    for batch in range(8):
        seed=int(np.random.SeedSequence([20261029209,tag,batch]).generate_state(1,dtype=np.uint64)[0]);seeds.append(seed);tick=time.perf_counter();pieces.append(k.stocks(d,pairs,p,x,configs,rates,seed));print(name,batch,'seconds',round(time.perf_counter()-tick,2),flush=True)
    stocks=combine(pieces);np.savez_compressed(folder/'stocks.npz',stocks=stocks,seeds=np.array(seeds,dtype=np.uint64),configs=configs);evaluate_stocks(stocks,np.load(LEVERAGE/'reference.npz')['stats'],arr(np.tile(x,(len(configs),1))),np.arange(len(configs),dtype=float),pd.read_csv(PARENT/'spx_membership.csv'),folder/'SPX');return relabel(folder/'SPX')

def run():
    OUT.mkdir(exist_ok=False);(OUT/'source').mkdir();sources=[Path(__file__)]+[HERE/n for n in ['m2_available_variance.cpp','m2_rich_protected_variance.cpp','m2_rich_protected_variance.py','m2_protected_variance.cpp','m2_protected_variance.py','asymmetric_readout.cpp','leverage_scan.cpp','two_speed_precision_kernel.cpp']]+[PARENT/'selection.json']
    for s in sources:shutil.copy2(s,OUT/'source'/s.name)
    dump(OUT/'source_hashes.json',{str(s):hashlib.sha256(s.read_bytes()).hexdigest() for s in sources});lib=OUT/'available.dylib';subprocess.run(['clang++','-O3','-std=c++17','-dynamiclib',str(HERE/'m2_available_variance.cpp'),'-o',str(lib)],check=True);prior,dated,_,na=context();assert np.allclose(prior['ends_years'],ENDSAV,atol=1e-14);x=arr(json.loads((PARENT/'selection.json').read_text())['parameters']);p=arr(prior['selected']);rates=arr(dated.layout.rates[MODES]);d=dict(np.load(HERE/'runs/m2_targeted_skew_001/grid_refinement/grid.npz'));k=Kernel(lib)
    configs=arr([[0,0,0,0],[0,0,0,1]]+[[j,l,eta*sign,2] for j in range(3) for l in range(2) for eta in [.25,.5,1.] for sign in [-1,1]])
    dump(OUT/'design.json',dict(configs=configs.tolist(),mode_indices=MODES.tolist(),half_life_days=(np.log(2)*365/rates).tolist(),baseline=x.tolist(),normalizer=na,calendar_ends_days=(ENDSAV*365).tolist(),damping_lambda=LAM,stress_only_availability=True,curve_refit=False,original_floor_retained_for_new_corrections=float((.05+.07*x[19])**2),no_OU_feedback=True,adopted=False))
    audit(k,d,p,x,rates);scores=stage(k,d,p,x,rates,configs,'pilot',32,111);best=min([s for s in scores if s['steps']==256 and s['candidate']>=2],key=lambda s:s['loss'])['candidate'];selected=arr([configs[0],configs[1],configs[best]]);dump(OUT/'selection.json',dict(pilot_new_best=best,configs=selected.tolist(),frozen_before_fresh=True,criterion='minimum all-primary pilot fine-grid SPX interval loss'))
    fresh=stage(k,d,p,x,rates,selected,'fresh',64,222);dump(OUT/'run.json',dict(status='pricing complete; scalar recurrence and reporting pending',pilot_paths=512,fresh_paths=1024,pilot_scores=scores,fresh_scores=fresh,adopted=False,production_baseline_unchanged=True));print('AVAILABLE PRICING COMPLETE',flush=True)

if __name__=='__main__':run()
