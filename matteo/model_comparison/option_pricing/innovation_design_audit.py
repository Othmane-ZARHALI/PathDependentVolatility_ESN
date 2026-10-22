"""Post-hoc design diagnostics only: covariance roundoff and old-pilot ablations."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .spx_innovation_study import datasets
from .spx_path_score_controls import ridge_fit
from esn_eq.config import ArchitectureConfig,RiskPremium
from esn_eq.pricing import DatedM2Pricer
from esn_eq.states import joint_transition_covariance,symmetric_square_root


def audit(output):
    record=json.loads((output/'run.json').read_text());quotes=pd.read_csv(output/'contracts.csv')
    p=DatedM2Pricer(ArchitectureConfig.from_mapping(record['architecture']),RiskPremium(**record['premium']))
    horizon=float(quotes.T_years.max());implemented={};covariance_rows=[]
    def represented(dt):
        if dt not in implemented:
            q=joint_transition_covariance(p.layout,float(dt));cross=q[1:,0]
            root=symmetric_square_root(q[1:,1:]-np.outer(cross,cross)/dt)
            implemented[dt]=root@root.T+np.outer(cross,cross)/dt
        return implemented[dt]
    with np.load(output/'normalizers.npz') as normalizers:
        fine=np.r_[normalizers['times3'],horizon]
        for level in (1,2,3):
            times=np.r_[normalizers[f'times{level}'],horizon];indices=np.searchsorted(fine,times)
            np.testing.assert_array_equal(fine[indices],times)
            patterns=set(tuple(np.diff(fine[a:b+1])) for a,b in zip(indices[:-1],indices[1:]))
            for pattern in patterns:
                accumulated=np.zeros((len(p.layout.rates),len(p.layout.rates)))
                for dt in pattern:
                    decay=np.exp(-p.layout.rates*dt);accumulated=decay[:,None]*accumulated*decay[None,:]+represented(dt)
                dt=float(sum(pattern));q=joint_transition_covariance(p.layout,dt)[1:,1:]
                exact=p.projection@q@p.projection.T;actual=p.projection@accumulated@p.projection.T
                covariance_rows.append(dict(level=level,fine_step_count=len(pattern),dt=dt,
                    max_projected_covariance_error=float(np.max(np.abs(actual-exact))),
                    max_relative_diagonal_error=float(np.max(np.abs(np.diag(actual-exact))/np.diag(exact)))))
    covariance=dict(description='Target projected covariance versus implemented fine conditional roots composed on each grid; numerical roundoff diagnostic, not a price bound.',rows=covariance_rows)
    (output/'covariance_roundoff_audit.json').write_text(json.dumps(covariance,indent=2)+'\n')
    dictionaries={'base':[], 'linear':[0,1,2], 'linear_square':[0,1,2,3,4,5], 'linear_mixed':[0,1,2,6,7,8], 'all9':list(range(9))}
    ablation=[]
    with np.load(output/'pilot_draws.npz') as pilot:
        data={key:pilot[key] for key in pilot.files}
        for target,(y,base,vol) in datasets(data).items():
            n=len(y);halves=(np.arange(n//2),np.arange(n//2,n))
            for fold in (0,1):
                training,validation=halves[fold],halves[1-fold];values={}
                for name,columns in dictionaries.items():
                    c=np.concatenate([base,vol[:,:,columns]],axis=2)
                    coefficients=np.stack([ridge_fit(y[training,j],c[training,j]) for j in range(y.shape[1])])
                    residual=(y[validation]-np.einsum('njk,jk->nj',c[validation],coefficients))*quotes.forward.to_numpy()/quotes.price_scale_U.to_numpy()
                    values[name]=float(residual[:,[2,3,8,9]].var(axis=0,ddof=1).mean())
                for name,value in values.items():
                    ablation.append(dict(target=target,fold=fold,dictionary=name,new_columns=len(dictionaries[name]),
                        central_scaled_variance=value,variance_gain_vs_base=values['base']/value))
    design=dict(status='post-hoc analysis of OLD training paths; no independent confirmation or change to reported checking estimates',
        purpose='Identify smaller dictionaries for a future freshly checked cost comparison; no subset runtime was measured here.',rows=ablation)
    (output/'pilot_ablation.json').write_text(json.dumps(design,indent=2)+'\n')
    print('Max composed covariance error:',max(row['max_projected_covariance_error'] for row in covariance_rows))
    print('Max relative diagonal error:',max(row['max_relative_diagonal_error'] for row in covariance_rows))
    for row in ablation:
        if row['target']=='price512':print(row)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path);audit(parser.parse_args().output)
