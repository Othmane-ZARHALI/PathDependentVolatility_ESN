from pathlib import Path
import json,hashlib,shutil,re
import numpy as np
import pandas as pd
h=Path('matteo/model_comparison/option_pricing');r=h/'runs/m2_stress_bank_002';parent=h/'runs/m2_stress_bank_001'
q=pd.read_csv(r/'fresh_check/prices_otm_bank.csv');groups=[g.index.to_numpy() for _,g in q.groupby('expiry_days',sort=True)]
def loss(p):
    dist=np.where(p<q.price_bid_U,p-q.price_bid_U,np.where(p>q.price_ask_U,p-q.price_ask_U,0))/q.price_scale_U
    return float(np.mean([np.mean(np.asarray(dist)[ids]**2) for ids in groups]))
with np.load(r/'fresh_check/blocks_bank.npz') as a,np.load(r/'fresh_check/blocks_two_speed.npz') as b:
    aa=a['native'];bb=b['native'];delta=loss(aa.mean(0))-loss(bb.mean(0));jk=[]
    for i in range(8):
        ids=np.arange(8)!=i;jk.append(loss(aa[ids].mean(0))-loss(bb[ids].mean(0)))
pair=dict(candidate='bank',reference='two_speed',option_loss_change=delta,paired_jackknife_se=float(np.std(jk,ddof=0)*np.sqrt(7)),leave_one_net_changes=jk)
(r/'fresh_check/paired_option_loss.json').write_text(json.dumps(pair,indent=2)+'\n')
f=h/'STRESS_BANK_RESULTS.md';s=f.read_text();s=s.replace('Three starts were attempted; all successes/failures', 'The winning candidate was the original converged third start, already matching forwards within8.5e-12; the finishing step did not change the selected model. Three starts were attempted; all successes/failures')
s=s.replace('## Adding further half-lives',f"Option-only paired checking loss change: {delta:+.6f}, jackknife SE {pair['paired_jackknife_se']:.6f}. The apparent small ranking difference should not be treated as established model superiority.\n\n## Adding further half-lives")
s=s.replace('- [Four smile panels]', '- [Native IV panels, using each model forward](visualizations/stress-bank025/vix-stress-bank-native-IV.png)\n- [Four smile panels]')
s=s.replace('canonical-call plots use the market forward.', 'canonical-call plots use the market forward. Because checking forward errors are material, their left-wing IVs can be unavailable even for a positive call estimate below market-forward intrinsic value. The separate native-IV panels use each model forward and avoid that parity-display issue; they still inherit the forward-estimation uncertainty.')
f.write_text(s)
for name in ['STRESS_BANK_RESULTS.md','plot_stress_bank_native.py']:shutil.copy2(h/name,r/'source'/name)
shutil.copy2(Path(__file__),r/'source/finalize_stress_bank.py')
for folder,lognames in [(parent,['m2-stress-bank-001.log','m2-stress-bank-001-audit.log']),(r,['m2-stress-bank-002.log','m2-stress-bank-audit.log'])]:
    logs=folder/'logs';logs.mkdir(exist_ok=True)
    for name in lognames:shutil.copy2(Path('/private/tmp')/name,logs/name)
    run=json.loads((folder/'run.json').read_text());run['status']='completed pilot; numerical audit PASS; forward population accuracy unresolved';run['audit']='audit.json';(folder/'run.json').write_text(json.dumps(run,indent=2)+'\n')
    manifest={'files':{str(p.relative_to(folder)):hashlib.sha256(p.read_bytes()).hexdigest() for p in folder.rglob('*') if p.is_file() and p.name!='artifact_manifest.json'},'figures':{str(p.relative_to(h)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (h/'visualizations/stress-bank025').glob('*')} if folder==r else {}}
    (folder/'artifact_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    for name,sha in manifest['files'].items():assert hashlib.sha256((folder/name).read_bytes()).hexdigest()==sha,name
    for name,sha in manifest['figures'].items():assert hashlib.sha256((h/name).read_bytes()).hexdigest()==sha,name
    print(folder.name,len(manifest['files']),'artifacts',len(manifest['figures']),'figures: hashes PASS')
for name in ['STRESS_BANK_PLAN.md','STRESS_BANK_RESULTS.md','STRESS_BANK_REFINEMENT.md']:
    s=(h/name).read_text();assert s.count('$$')%2==0,name
    for link in re.findall(r'\]\(([^)]+)\)',s):
        if '://' not in link:assert (h/link).exists(),(name,link)
print('Markdown math and local links PASS');print('Paired option loss:',json.dumps(pair))
