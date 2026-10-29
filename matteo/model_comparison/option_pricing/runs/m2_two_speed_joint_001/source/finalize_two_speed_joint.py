from pathlib import Path
import json,hashlib,re,shutil
import numpy as np
import pandas as pd
h=Path('matteo/model_comparison/option_pricing');r=h/'runs/m2_two_speed_joint_001'
# Preserve additional raw weighted-tail diagnostics; never substitute prices.
rows=[]
for name in ['vix_only','joint']:
 q=pd.read_csv(r/'checking'/f'prices_spx_{name}.csv')
 with np.load(r/'checking'/f'stocks_{name}.npz') as a:
  for day,g in q.groupby('expiry_days',sort=True):
   S=a[f'd{day}_stock'];w=a[f'd{day}_weights'];values=a[f'd{day}_values']
   for j,z in enumerate(g.itertuples()):
    sign=1 if z.option_type=='call' else -1;pay=w*np.maximum(sign*(S-z.K_over_F),0)*z.forward;total=pay.sum();top=np.sort(pay)[::-1];sq=(values[:,j]-values[:,j].mean())**2
    rows.append(dict(candidate=name,day=int(day),strike=float(z.strike),option_type=z.option_type,positive_raw_payoffs=int((pay>0).sum()),weighted_hit_probability=float(np.mean(w*(pay>0))),raw_price=float(pay.mean()),payoff_ESS=float(total**2/(pay@pay)) if total else 0.,top10_payoff_fraction=float(top[:10].sum()/total) if total else 0.,top10_controlled_centered_square_fraction=float(np.sort(sq)[-10:].sum()/sq.sum()) if sq.sum() else 0.))
pd.DataFrame(rows).to_csv(r/'SPX_weighted_tail_diagnostics.csv',index=False)
# Remove a dead plotting expression and a harmless generator escape warning.
f=h/'plot_two_speed_joint.py';s=f.read_text();s=s.replace(";p=np.array(json.loads((run/'source/run.json').read_text())['selected']) if False else None",'');f.write_text(s)
generator=Path('/private/tmp/report_two_speed_joint.py');s=generator.read_text().replace("'**SPX accuracy gates fail:**", "r'**SPX accuracy gates fail:**");generator.write_text(s)
f=h/'TWO_SPEED_JOINT_RESULTS.md';s=f.read_text().replace('See control_variance_diagnostic.csv.', 'See control_variance_diagnostic.csv. Raw likelihood-weighted payoff hits, ESS, concentration and controlled-estimator variance concentration for every primary SPX quote are retained in SPX_weighted_tail_diagnostics.csv; zero sampled payoffs do not certify a zero theoretical price.')
f.write_text(s)
for f in [h/'plot_two_speed_joint.py',h/'TWO_SPEED_JOINT_RESULTS.md',generator,Path(__file__)]:shutil.copy2(f,r/'source'/f.name)
# Record derived-cache integrity separately from the original-input audit.
c=h/'runs/m2_two_speed_joint_cache_001';cm={str(p.relative_to(c)):hashlib.sha256(p.read_bytes()).hexdigest() for p in c.rglob('*') if p.is_file() and p.name!='artifact_manifest.json'};(c/'artifact_manifest.json').write_text(json.dumps(dict(files=cm),indent=2)+'\n')
extra=h/'runs/m2_maturity_curve_fit_001/convolution.dylib';shutil.copy2(extra,r/'source/frozen_normalizer_convolution.dylib')
(r/'auxiliary_input_hashes.json').write_text(json.dumps({str(extra):hashlib.sha256(extra.read_bytes()).hexdigest()},indent=2)+'\n')
manifest={'files':{str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in r.rglob('*') if p.is_file() and p.name!='artifact_manifest.json'},'figures':{str(p.relative_to(h)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (h/'visualizations/two-speed-joint026').glob('*')}};(r/'artifact_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
for name,sha in manifest['files'].items():assert hashlib.sha256((r/name).read_bytes()).hexdigest()==sha,name
for name,sha in manifest['figures'].items():assert hashlib.sha256((h/name).read_bytes()).hexdigest()==sha,name
for name in ['TWO_SPEED_JOINT_RESULTS.md','TWO_SPEED_JOINT_PLAN.md']:
 s=(h/name).read_text();assert s.count('$$')%2==0,name
 for target in re.findall(r'\]\(([^)]+)\)',s):
  if '://' not in target:assert (h/target).exists(),(name,target)
print('VERIFIED',len(manifest['files']),'run artifacts',len(manifest['figures']),'figures;',len(cm),'derived cache files; Markdown math/links PASS')
