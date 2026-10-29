from pathlib import Path
import re,hashlib,json,zipfile,subprocess,sys
root=Path.cwd();d=root/'docs/research';m=json.loads((d/'EVIDENCE_MANIFEST.json').read_text());f=json.loads((d/'FIGURE_PROVENANCE.json').read_text());problems=[]
def sha(b):return hashlib.sha256(b).hexdigest()
for e in m['reports']:
 p=root/e['path']
 if p.name in ['CONTINUE_FITTING_PROMPT.md','NEXT_STEPS.md'] and p.parent.name=='option_pricing':
  original=(d/'archive/2026-10-29'/p.name).read_bytes()
 else:
  original=p.read_bytes()
  if original.startswith(b'<!-- Research consolidation 2026-10-29 -->'):
   original=original.split(b'\n\n',1)[1]
 if sha(original)!=e['sha256']:problems.append('Original report hash '+e['path'])
for e in f['figures']:
 if sha((root/e['source']).read_bytes())!=e['sha256'] or sha((root/e['copy']).read_bytes())!=e['sha256']:problems.append('Figure changed '+e['source'])
for e in m['visual_assets']:
 if sha((root/e['path']).read_bytes())!=e['sha256']:problems.append('Visual source changed '+e['path'])
if sha((root/m['scope']['path']).read_bytes())!=m['scope']['sha256']:problems.append('Scope changed')
# Check links in new active prose and the two replacement entry points.
files=list(d.glob('*.md'))+[root/'matteo/model_comparison/option_pricing'/n for n in ['NEXT_STEPS.md','CONTINUE_FITTING_PROMPT.md']]
link_count=0
for p in files:
 for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
  if '://' in target:continue
  name,_,anchor=target.partition('#');q=(p.parent/name).resolve() if name else p
  link_count+=1
  if not q.exists():problems.append('Broken link '+str(p.relative_to(root))+': '+target)
  elif anchor and q.suffix=='.md':
   headings=[]
   for title in re.findall(r'^#+\s+(.+)$',q.read_text(),re.M):
    headings.append(re.sub(r'[^\w\- ]','',title.lower()).replace(' ','-'))
   if anchor not in headings:problems.append('Broken anchor '+target)
recovery=json.loads((d/'archive/2026-10-29/RECOVERY_MANIFEST.json').read_text())
recovery_zip=d/'archive/2026-10-29/missing_tracked_evidence_HEAD.zip'
if recovery_zip.exists() and '--git-history' not in sys.argv:
 with zipfile.ZipFile(recovery_zip) as z:
  for e in recovery['files']:
   if sha(z.read(e['path']))!=e['sha256']:problems.append('Recovered hash '+e['path'])
else:
 # These committed ancestor bytes need not be uploaded again as a large ZIP.
 proc=subprocess.Popen(['git','cat-file','--batch'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
 try:
  for e in recovery['files']:
   proc.stdin.write((recovery['git_commit']+':'+e['path']+'\n').encode());proc.stdin.flush()
   header=proc.stdout.readline().split()
   if len(header)!=3 or header[1]!=b'blob':
    problems.append('Historical Git object missing '+e['path']);continue
   size=int(header[2]);parts=[];remaining=size
   while remaining:
    chunk=proc.stdout.read(min(remaining,1048576))
    if not chunk:raise RuntimeError('Truncated historical Git blob')
    parts.append(chunk);remaining-=len(chunk)
   proc.stdout.read(1)
   if sha(b''.join(parts))!=e['sha256']:problems.append('Historical Git hash '+e['path'])
 finally:
  proc.stdin.close();proc.stdout.close();proc.wait()
# Verify pre-existing deletions have not been restored or otherwise changed.
old_deleted=[e['path'] for e in recovery['files']]
if any((root/s).exists() for s in old_deleted):problems.append('Pre-existing deletion restored')
result={'reports_preserved':len(m['reports']),'selected_figures_unchanged':len(f['figures']),'visual_source_assets_unchanged':len(m['visual_assets']),'local_links_checked':link_count,'recovered_HEAD_files_verified':len(recovery['files']),'missing_referenced_run_directories':len(m['referenced_run_directories_missing']),'referenced_present_runs_without_run_json':len(m['referenced_present_runs_without_run_json']),'scope_unchanged':True,'problems':problems}
(d/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if problems:raise SystemExit(1)
