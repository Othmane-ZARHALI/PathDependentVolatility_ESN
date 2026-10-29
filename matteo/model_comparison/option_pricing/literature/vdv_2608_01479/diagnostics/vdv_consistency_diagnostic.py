"""Independent algebraic checks for a VIX-first inverse problem; no market fit."""
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.integrate import quad_vec
from scipy.linalg import expm
from .full_readout_fit import HERE,dump
OUT=HERE/'literature/vdv_2608_01479/diagnostics'

def run():
    OUT.mkdir(exist_ok=False);D=30/365
    # Coefficients transcribed and visually verified from Section4.3, PDFpage19.
    displayed=np.array([-101.37,-4.19,1.16]);roots=np.sort(np.roots(displayed[::-1]));positive=float(roots[roots>0][0]);assert np.polynomial.polynomial.polyval(10,displayed)<0;assert np.polynomial.polynomial.polyval(20,displayed)>0
    # Our own illustrative parameters, not the authors' calibrated local-volatility model.
    k,theta,sigma=12.,.2,1.;L=np.array([[0,k*theta,0],[0,-k,2*k*theta],[0,0,sigma*sigma-2*k]])
    B=np.zeros((6,6));B[:3,:3]=L;B[:3,3:]=np.eye(3);A=expm(B*D)[:3,3:]/D;c=np.linalg.solve(A,np.array([0.,0.,1.]));A2,error=quad_vec(lambda u:expm(L*u)/D,0,D,epsabs=1e-12,epsrel=1e-12);matrix_error=float(abs(A-A2).max());coefficient_error=float(abs(A@c-[0,0,1]).max());probe=np.linspace(.005,.5,1001);rolling=np.column_stack([np.ones(len(probe)),probe,probe**2])@(A2@c);residual=float(abs(rolling-probe**2).max());xstar=float(-c[1]/(2*c[2]));minimum=float(np.polynomial.polynomial.polyval(xstar,c));assert matrix_error<1e-12 and coefficient_error<1e-12 and residual<1e-12;assert minimum<0
    # Closed-form positive control with autonomous geometric-Brownian VIX.
    mu,vol=.1,1.;lam=2*mu+vol*vol;factor=float(D*lam/np.expm1(D*lam));oracle=quad_vec(lambda u:factor*np.exp(lam*u)/D,0,D)[0];assert abs(oracle-1)<1e-12
    result=dict(status='completed independent algebraic diagnostics; not calibration',source=dict(url='https://arxiv.org/abs/2608.01479',version='v1',paper_pdf_sha256=hashlib.sha256((OUT.parent/'paper.pdf').read_bytes()).hexdigest(),polynomial_source='Section4.3, PDFpage19; visually checked; rounded coefficients'),displayed_terminal_polynomial=dict(coefficients_ascending=displayed.tolist(),real_roots=roots.tolist(),negative_region_on_positive_axis=[0,positive],p_at10=float(np.polynomial.polynomial.polyval(10,displayed)),units='x in the displayed numerical example uses VIX points; positive scaling of variance does not alter its sign',scope='sign check of displayed terminal polynomial only; actual support, boundary extrapolation and use of a floor in author code not established'),illustrative_Brennan_Schwartz=dict(parameters=dict(kappa=k,theta=theta,sigma=sigma,window=D),generator_on_coefficients=L.tolist(),rolling_operator=A.tolist(),coupling_coefficients_ascending=c.tolist(),negative_interval=np.sort(np.roots(c[::-1])).tolist(),minimum=minimum,minimizer=xstar,block_exponential_vs_independent_quadrature_error=matrix_error,coefficient_residual=coefficient_error,grid_rolling_identity_error=residual,scope='a polynomial satisfies the rolling identity algebraically but is not a globally nonnegative variance; no impossibility claim for all time-dependent couplings'),positive_GBM_control=dict(mu=mu,sigma=vol,window=D,coupling_coefficient=factor,rolling_identity_error=float(abs(oracle-1))),floor_identity='If psi is exactly consistent, flooring to max(psi,0) adds E_t[integral(-psi)_+]/Delta to implied squared VIX. This is positive when occupation of the negative region has positive conditional probability and magnitude.',no_floor_applied=True,baseline_changed=False,no_new_Lean_claim=True)
    dump(OUT/'algebraic_audit.json',result);print('VDV ALGEBRAIC AUDIT PASS',positive,matrix_error,residual,minimum,flush=True)
if __name__=='__main__':run()
