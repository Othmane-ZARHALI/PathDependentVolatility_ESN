"""VIX-focused curve fits; two full SPX trade-off evaluations, no adoption."""
import argparse,ctypes as ct,json,time
from pathlib import Path
import numpy as np,pandas as pd
from scipy.optimize import least_squares
from .cap_tail_study import engine,tails,VARIANTS
from .full_readout_fit import HERE,dump
from .m2_spx import ROOT,frozen_stage3
from .maturity_curve import CurveKernel
from .fit_parity_parameters_v3 import save_result
from .provenance import capture_sources,save_sources,source_hashes,check_sources
from esn_eq.pricing import DatedM2Pricer


def focused_fit(e,p):
    # z_i=(xi_i-f)/(.2-f); exact hard constraints xi>=f and xi<=.2.
    z0=np.r_[(p[7:-1]-p[-1])/(.2-p[-1]),(p[-1]-.0025)/.0011];trials=[]
    def unpack(z):
        x=p.copy();x[-1]=.0025+.0011*z[-1];x[7:-1]=x[-1]+(.2-x[-1])*z[:-1];return x
    def calc(z):
        x=unpack(z);v=e.vix_prediction(x);native=v['prices'].copy();index=0
        for j,(_,g) in enumerate(e.vix.groupby('expiry_days',sort=True)):
            n=len(g);native[index:index+n]-=g.option_type.eq('put').to_numpy()*(v['forwards'][j]-g.strike.to_numpy());index+=n
        r=e.otm_loss.evaluate(native,v['forwards'],x[7:-1]);trials.append(dict(parameters=x.tolist(),loss=float(r['total']),components=r['components']))
        return r['residuals']
    fit=least_squares(calc,z0,bounds=(np.zeros(17),np.ones(17)),max_nfev=40,ftol=1e-8,xtol=1e-8,gtol=1e-8)
    return unpack(fit.x),dict(success=bool(fit.success),termination=fit.message,nfev=int(fit.nfev),actual_residual_calls=len(trials),trials=trials,method='VIX-only16-curve-plus-floor fit, readout fixed; no joint safeguards; same floor5–6% comparison range; no global optimum claim')


def run(out):
    sub=out/'followup';sub.mkdir(exist_ok=False);prior=HERE/'runs/m2_maturity_curve_fit_001';d=json.loads((prior/'run.json').read_text());p=np.array(d['selected']);ends=np.array(d['ends_years']);cache=Path(d['cache']);spx=pd.read_csv(cache/'spx_contracts.csv');vix=pd.read_csv(cache/'vix_contracts.csv');f=pd.read_csv(cache/'forward_targets.csv');vc=json.loads((prior/'vix_controls.json').read_text())
    sources=[Path(__file__),HERE/'cap_tail_study.py',HERE/'cap_tail_kernel.cpp',out/'scores.json',out/'design.json',out/'cap.dylib'];captured=capture_sources(sources);save_sources(sub,captured,ROOT)
    dump(sub/'design.json',dict(selection=['higher_L','logtail025'],reason='TRAINING profiles show steeper right wings; higherL tests strong slope/forward tradeoff; logtail025 tests nonsaturating growth with nearly unchanged native loss. No checking prices used.',max_fits=2,max_nfev_per_fit=40,full_stock_evaluations=2,source_sha256=source_hashes(captured,ROOT)))
    k=CurveKernel(out/'cap.dylib',out/'cap.dylib',ends);k.lib.configure.argtypes=[ct.c_int,ct.c_double];a,_,premium,_=frozen_stage3();dated=DatedM2Pricer(a,premium);initial=np.array(d['initial_state']);scores=[];tailrows=[];mom=[]
    for name in ['higher_L','logtail025']:
        _,mode,strength,dl,kap=next(x for x in VARIANTS if x[0]==name);k.lib.configure(mode,strength);x=p.copy();x[5]+=dl
        e=engine(out,cache,k,dated,initial,spx,vix,f,vc,Path(d['controls']),d)
        print('FOCUSED FIT',name,flush=True);tick=time.perf_counter();selected,fit=focused_fit(e,x);dump(sub/f'fit_{name}.json',fit);dump(sub/f'parameters_{name}.json',dict(profile=x.tolist(),fitted=selected.tolist()));e.close()
        # Freeze before separate checking and actual SPX tradeoff prices.
        print('FULL STOCK',name,time.perf_counter()-tick,flush=True);full=engine(out,cache,k,dated,initial,spx,vix,f,vc,Path(d['controls']),d,False)
        try:
            r=full.predict(selected,True,False);scores.append(save_result(full,r,selected,name+'_fit',sub));print('FULL SCORE',name,r['components'],r['original_otm_components'],flush=True)
        finally:full.close()
        cc=Path(d['check_cache']);ce=engine(out,cc,k,dated,initial,spx,vix,f,vc,Path(d['check_controls']),d)
        try:
            for label,x0 in [(name,x),(name+'_fit',selected)]:
                print('CHECK VIX',label,flush=True);rr=ce.predict(x0,True,False);folder=sub/'checking_vix';folder.mkdir(exist_ok=True);save_result(ce,rr,x0,label,folder);t,m=tails(folder,label,vix,x0,np.array(d['ends_days']));tailrows+=t;mom+=m
        finally:ce.close()
        dump(sub/'scores.json',scores);pd.DataFrame(tailrows).to_csv(sub/'checking_tails.csv',index=False);pd.DataFrame(mom).to_csv(sub/'checking_moments.csv',index=False)
    check_sources(captured);dump(sub/'run.json',dict(status='completed two VIX-focused fits, all-seven-maturity SPX tradeoffs and reused separate VIX checking',accepted_joint_calibrations=0,scores=scores))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);run(p.parse_args().run.resolve())
