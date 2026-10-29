"""Independent scalar quadrature, alpha-zero replay, and paired score diagnostics."""
import json,time
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.special import expit
from .slow_readout import *
from .asymmetric_readout import Kernel as OriginalKernel
from .two_speed_joint_core import physical

def run():
    prior,dated,norm,na=context();p=arr(prior['selected']);initial=arr(prior['initial_state']);design=json.loads((OUT/'design.json').read_text());k=Kernel(OUT/'slow.dylib');errors=[];correlations=[]
    for idx in design['mode_indices']:
        law=Law(dated,idx);m,c=law(27/365,initial);m=m+arr([.1,-.2,.3,-.4]);alpha=arr([0,.5,1]);order=24;z,w=hermgauss(order);z*=np.sqrt(2);w/=np.sqrt(np.pi);yy,ff=np.meshgrid(z,z,indexing='ij');root=symmetric_square_root(c[1:3,1:3]);off=np.column_stack([yy.ravel(),ff.ravel()])@root.T;weights=np.outer(w,w).ravel();beta=c[3,1:3]@np.linalg.pinv(c[1:3,1:3],rcond=1e-12);sd=np.sqrt(max(0,c[3,3]-beta@c[1:3,3]));Y=m[1]+off[:,0];F=m[2]+off[:,1];eta=p[0]*(m[0]+np.sqrt(c[0,0])*z[None,:])+p[1]*(Y[:,None]-p[2])**2+p[3]*(F[:,None]-p[4])**2;G=np.exp(2*(p[5]-np.logaddexp(0,p[6]*(p[5]-eta))/p[6]));Eg=G@w;center=float(law(27/365,initial)[0][3]);Xm=m[3]+off@beta-center;D=np.exp(np.logaddexp(0,-(Xm[:,None]+sd*z[None,:]))[:,:,None]*alpha[None,None,:]);ED=(D*w[None,:,None]).sum(1);oracle=(weights[:,None]*Eg[:,None]*ED).sum(0);fast=k.moments4(m,c,p,center,alpha,order)[0,:,1];errors.append(float(abs(fast/oracle-1).max()));correlations.append((c[3,1:3]/np.sqrt(c[3,3]*np.diag(c)[1:3])).tolist())
    # Independent existing implementation must exactly recover alpha-zero stock paths.
    d=dict(np.load(HERE/'runs/m2_targeted_skew_001/grid_refinement/grid.npz'));x=arr(design['baseline']);models=arr([x]);new=k.stocks(d,1,p,models,arr([0]),np.ones((len(d['dt']),1,2)),[design['mode_indices'][0]],20261028614);oldk=OriginalKernel(OUT/'slow.dylib');ref=oldk.stocks(d,1,p,models,arr([0]),np.ones((len(d['dt']),1)),1.,20261028614);replay=float(abs(new-ref).max());assert replay<1e-11
    assert max(errors)<1e-5
    dump(OUT/'implementation_audit.json',dict(scalar_conditional_tensor_GH24_relative_errors=errors,slow_correlations_with_Y_F=correlations,alpha_zero_existing_stock_max_error=replay,normalization_identity='E_old[(aR+wA)D*(a+wA)/(a*mRD+wA*mD)]=a+wA',mathematical_scope='conditional identity given independent stress states; no full VIX-law invariance',code_scope='numerical checks, no Lean certification',all_quotes_retained=True,old_fast_readout_retained=True))
    print('AUDIT',errors,replay,flush=True)
if __name__=='__main__':run()
