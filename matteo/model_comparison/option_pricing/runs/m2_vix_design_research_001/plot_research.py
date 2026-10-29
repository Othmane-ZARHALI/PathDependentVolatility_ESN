"""Static evidence figures; no new model simulations or digitized Othmane curves."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matteo.model_comparison.option_pricing.plot_maturity_curve import iv
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
cert=json.loads((HERE/'moment_certificates.json').read_text());rows=[r for r in cert['records'] if r['forward_tolerance']==0]
fig,ax=plt.subplots(figsize=(10,6),layout='constrained');x=np.arange(4);w=.35;current=[float(__import__('fractions').Fraction(r['baseline_second_target'])) for r in rows];required=[r['certified_lower_bound_float'] for r in rows]
a=ax.bar(x-w/2,current,w,label='Frozen curve budget',color='#2865ad');b=ax.bar(x+w/2,required,w,label='Certified necessary lower bound from VIX quotes',color='#d66b36')
ax.bar_label(a,fmt='%.1f',padding=3);ax.bar_label(b,fmt='%.1f',padding=3);ax.set(xticks=x,xticklabels=[str(r['day']) for r in rows],xlabel='VIX expiry (days)',ylabel='E[VIX²] in VIX points²',ylim=(0,790),title='The frozen variance budget is insufficient at three VIX expiries\nAll quoted bid–ask intervals and observed VIX forwards');ax.legend(frameon=False,loc='upper left');ax.grid(axis='y',alpha=.18)
for ext in ['png','svg']:fig.savefig(HERE/f'variance-budget.{ext}',dpi=150,bbox_inches='tight')
plt.close(fig)
prior=ROOT/'runs/m2_postcap_vov_001/fresh_check';tables={name:pd.read_csv(prior/f'prices_canonical_{name}.csv') for name in ['baseline','g1_h7']}
fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained');ranges={27:(-.26,.91),55:(-.29,.98),118:(-.28,.95),153:(-.25,.93)}
for ax,day in zip(axes.flat,[27,55,118,153]):
 q=tables['baseline'];q=q[q.expiry_days.eq(day)].sort_values('strike');lo,hi=ranges[day]
 ax.axvspan(np.exp(lo),np.exp(hi),color='#c9c9c9',alpha=.3,label='cs_doc calibration strike range');ax.plot(q.K_over_F,iv(q,q.price_mid_U),'k-',label='Market, full retained quotes')
 for name,color,label in [('baseline','#2865ad','Our original baseline'),('g1_h7','#168958','Our tested post-cap multiplier')]:
  g=tables[name];g=g[g.expiry_days.eq(day)].sort_values('strike');ax.plot(g.K_over_F,iv(g,g.model_price_U),color=color,label=label)
 ax.set(xscale='log',title=f'{day} days',xlabel='Strike / market forward',ylabel='Black call IV (%)');ax.grid(alpha=.18)
h,l=axes.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncols=2,frameon=False);fig.suptitle('Othmane comparison: the fitted strike range is materially narrower\nShading from cs_doc Table1; displayed model curves are our archived simulations')
for ext in ['png','svg']:fig.savefig(HERE/f'strike-range-comparison.{ext}',dpi=150,bbox_inches='tight')
plt.close(fig)
