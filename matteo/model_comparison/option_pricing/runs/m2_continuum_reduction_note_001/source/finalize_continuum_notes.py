from pathlib import Path
import hashlib,json,re,shutil
b=Path('/Users/matteogambara/Desktop/code/PathDependentVolatility_ESN/matteo/model_comparison/option_pricing');out=b/'runs/m2_continuum_reduction_note_001';frozen=b/'runs/m2_continuum_response_003';notes=[b/'CONTINUUM_EMPIRICAL_REDUCTION.md',b/'VIX_27_DAY_DIAGNOSIS.md']
# Check derived statements against saved numbers and current immutable experiment.
d=json.loads((out/'diagnostics.json').read_text());score=json.loads((frozen/'fresh_check/score_continuous.json').read_text());metrics=[r for r in d['shape_metrics'] if r['candidate']=='continuous'];assert abs(sum(r['native_option_loss_contribution'] for r in metrics)-score['option_loss'])<1e-12
assert d['market_snapshot_date_inferred_from_expiry_dates']=='2026-02-19';assert abs(d['second_moment_budget']['model_second_moment']-467.4495923035969)<1e-9
for p in notes:
 s=p.read_text();assert s.count('$$')%2==0,p;assert s.replace('$$','').count('$')%2==0,p;assert s.count('{')==s.count('}'),p
 for link in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',s):
  if not link.startswith(('https://','http://','#')):assert (p.parent/link.split('#')[0]).exists(),(p,link)
# Preserve and verify the archived calibration and all of its figures.
m=json.loads((frozen/'artifact_manifest.json').read_text())
for root,key in [(frozen,'sha256'),(b,'figure_sha256')]:
 for rel,h in m[key].items():assert hashlib.sha256((root/rel).read_bytes()).hexdigest()==h,rel
checks=dict(status='PASS',note_math_delimiters_and_local_links=True,derived_loss_sum=True,archived_calibration_hashes_unchanged=len(m['sha256']),archived_figures_unchanged=len(m['figure_sha256']),no_calibration_or_model_change=True)
(out/'documentation_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
checkpoint='''## 27 October 2026 — teaching notes and first VIX maturity diagnosis

Saved CONTINUUM_EMPIRICAL_REDUCTION.md: equations, Bergomi attribution, the Beta-to-two-atom limit, bounded empirical reduction and why the general integral remains useful for future market snapshots. The reduction concerns mixture/readout complexity, not stochastic state dimension; only one stored market snapshot was tested.

Saved VIX_27_DAY_DIAGNOSIS.md with reproducible read-only diagnostics in runs/m2_continuum_reduction_note_001. The27day right wing is too low while central calls are too high. In the matching boundary reference, the strike100 payoff-weighted stress variance share is98.71%; stress dominance and shared response parameters are a leading shape mechanism. Exact-forward necessary second moment468.99933 exceeds selected467.44959. A27day-only constrained calibration is a proposed discriminating experiment, not an authorized scope change or a completed fit. All217quotes/fourforwards, the working candidate and SPX pause are unchanged.

'''
for name in ['RESEARCH_LEDGER.md','NEXT_STEPS.md','CONTINUE_FITTING_PROMPT.md']:
 p=b/name
 if p.exists() and 'teaching notes and first VIX maturity diagnosis' not in p.read_text():p.write_text(checkpoint+p.read_text())
p=b/'CONTINUUM_RESPONSE_RESULTS.md';s=p.read_text();addition='\n## Subsequent explanatory notes\n\nSee [the empirical reduction and its mathematical meaning](CONTINUUM_EMPIRICAL_REDUCTION.md) and [the27day wing diagnosis](VIX_27_DAY_DIAGNOSIS.md). These add documentation and read-only diagnostics; the archived calibration, model parameters and original audit manifests are unchanged.\n';
if '## Subsequent explanatory notes' not in s:p.write_text(s+addition)
post=out/'source';post.mkdir(exist_ok=True)
for p in [*notes,b/'diagnose_vix27_continuum.py',Path('/private/tmp/finalize_continuum_notes.py')]:shutil.copy2(p,post/p.name)
manifest={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file() and p.name!='artifact_manifest.json'}
(out/'artifact_manifest.json').write_text(json.dumps(dict(sha256=manifest,count=len(manifest),scope='New documentation and read-only diagnostics; archived calibration preserved'),indent=2)+'\n')
print(json.dumps(checks));print('Saved',len(manifest),'new artifact hashes');print('Note words:',{p.name:len(p.read_text().split()) for p in notes})
