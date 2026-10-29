"""Deterministic curve/window illustration; no option prices or new fit."""
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def render():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    state = HERE/'runs/m2_joint_state_fit_003/selected.json'
    curve = HERE/'runs/m2_curve_refinement_003/selected.json'
    design = HERE/'runs/m2_curve_refinement_003/design.json'
    s, c, d = [json.loads(p.read_text()) for p in [state, curve, design]]
    out = HERE/'visualizations/variance-pieces-006'
    out.mkdir(parents=True, exist_ok=True)
    ends5 = np.array(d['ends_years'])*365
    ends3 = np.array([.25,.75,d['ends_years'][-1]])*365
    levels = [np.array(s['levels']), np.array(c['levels'])]
    windows = np.array([27,55,118,153])
    starts5 = np.r_[0,ends5[:-1]]
    overlap = np.maximum(0,np.minimum(windows[:,None]+30,ends5)-np.maximum(windows[:,None],starts5))/30
    fig, axes = plt.subplots(2,1,figsize=(12,8),gridspec_kw={'height_ratios':[2.1,1.3]},layout='constrained')
    for xi,ends,color,label,style in zip(levels,[ends3,ends5],['#b66d16','#246cb1'],
            ['Selected 3-piece curve','Selected 5-piece curve'],['--','-']):
        x=np.r_[0,ends[ends<365],365]
        y=np.r_[100*np.sqrt(xi),100*np.sqrt(xi[-1])]
        axes[0].step(x,y,where='post',lw=2,color=color,ls=style,label=label)
    for day in ends5[:-1]:axes[0].axvline(day,color='0.7',lw=.7)
    for left,right,xi in zip(starts5,ends5,levels[1]):
        midpoint=(left+min(right,365))/2
        axes[0].text(midpoint,100*np.sqrt(xi)+.7,f'{100*np.sqrt(xi):.2f}%',ha='center',fontsize=10,color='#246cb1')
    axes[0].set(xlim=(0,365),ylim=(15,31),xlabel='Future calendar day',ylabel='100 × square root of annual variance',
                title='Extra levels let different future periods move separately')
    axes[0].legend(frameon=False,loc='upper right');axes[0].grid(axis='y',alpha=.2)
    for j,day in enumerate(windows):
        axes[1].broken_barh([(day,30)],(j-.2,.4),facecolors='#468b81')
        axes[1].text(day+15,j+.27,f'{day}–{day+30}',ha='center',fontsize=10)
    for day in ends5[:-1]:axes[1].axvline(day,color='0.6',lw=.8)
    axes[1].set(xlim=(0,365),ylim=(-.6,3.7),xlabel='Future calendar day',yticks=range(4),
                yticklabels=[f'VIX expiry: {day} days' for day in windows],title='Each VIX expiry averages the following 30 days of variance')
    axes[1].invert_yaxis();axes[1].grid(axis='x',alpha=.15)
    fig.suptitle('Why test five pieces?\nSearching-selected curves; extra five-piece gain is small and inconclusive on fresh checks.',fontsize=14)
    fig.supxlabel('The percentages above are variance-input scales, not option implied volatilities. The final level continues beyond the first year.',fontsize=10)
    outputs=[]
    for ext in ['png','svg']:
        p=out/f'variance-pieces-and-vix-windows.{ext}';fig.savefig(p,dpi=150,bbox_inches='tight');outputs.append(p)
    plt.close(fig)
    (out/'window_weights.json').write_text(json.dumps(dict(expiries_days=windows.tolist(),ends_days=ends5.tolist(),
        fraction_of_each_window=overlap.tolist(),meaning='Population mean-square VIX identity weights for deterministic curve, not conditional path exposures or VIX forward predictions.'),indent=2)+'\n')
    (out/'manifest.json').write_text(json.dumps(dict(source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),state,curve,design]},
        output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs},new_prices=False,new_fit=False),indent=2)+'\n')
    (out/'README.md').write_text('''# Reading the curve picture

The upper panel compares the searching-selected three- and five-piece variance
inputs at the same initial state. Each percentage is the square root of an annual
variance level. It is not an option implied volatility or average volatility.
Fresh checks find only a 0.64% pooled extra gain with one trial reversal and
uncertainty larger than that gain. Three pieces remain the working reference;
neither is a qualified market fit.

The lower panel shows the next 30 days used by each VIX expiry. Splitting near
day 57 separates most of the first two windows; splitting near day 151 separates
the two later windows. Existing boundaries at days 91.25 and 273.75 remain.
The actual added second boundary is day 151.000543, on the unchanged stock grid.

SPX options and VIX options/forwards jointly constrain these levels. The last
piece starts at day 273.75 and is weakly constrained by the current fitting
expiries: only the 302-day SPX slice reaches it. Other retained SPX maturities
keep their reserved/diagnostic roles and were not added to the objective.
More pieces can adjust the timing of average variance, but cannot freely adjust
every strike's smile shape. No prices are changed or calculated in this picture.
''')


if __name__ == '__main__':
    render()
