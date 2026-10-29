"""Inspect model-native versus common-market-forward IV, without repricing."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def render(directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    q=pd.read_csv(directory/'iv_conventions.csv')
    names=list(q.candidate.unique())
    assert len(names)==2 and set(q.expiry_days)=={27,55,118,153}
    labels={names[0]:'Reference',names[1]:'Selected candidate'}
    colors={names[0]:'#c17621',names[1]:'#246cb1'}
    styles={names[0]:'--',names[1]:'-'}
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(4,2,figsize=(13,13),layout='constrained')
    for row,days in enumerate([27,55,118,153]):
        subset=q[q.expiry_days==days]
        allvalues=np.r_[subset.market_bid_iv,subset.market_ask_iv,subset.model_iv_market_forward,subset.model_iv_own_forward]
        lo,hi=np.nanmin(allvalues),np.nanmax(allvalues);pad=.08*(hi-lo)
        for col,convention in enumerate(['market','own']):
            ax=axes[row,col];market=subset[subset.candidate==names[0]].sort_values('strike')
            for side in ['put','call']:
                m=market[market.option_type==side]
                ax.fill_between(m.K_over_market_F,m.market_bid_iv,m.market_ask_iv,color='0.3',alpha=.15,label='Market bid–ask' if side=='put' else None)
                ax.plot(m.K_over_market_F,m.market_mid_iv,color='0.2',lw=1)
                for name in names:
                    g=subset[(subset.candidate==name)&(subset.option_type==side)].sort_values('strike')
                    field='model_iv_market_forward' if convention=='market' else 'model_iv_own_forward'
                    ax.plot(g.K_over_market_F,g[field],styles[name],color=colors[name],lw=1.7,label=labels[name] if side=='put' else None)
            ax.set_xscale('log');xs=market.K_over_market_F.to_numpy();ticks=np.linspace(xs.min(),xs.max(),4)
            ax.set_xticks(ticks,[f'{v:.2f}' for v in ticks]);ax.set_xticks([],minor=True)
            ax.set_ylim(lo-pad,hi+pad);ax.yaxis.set_major_locator(MaxNLocator(5))
            ax.axvline(1,color='0.6',lw=.7);ax.grid(alpha=.15)
            ax.set_xlabel('Strike / MARKET forward (log scale)');ax.set_ylabel('Annualized Black IV (%)')
            title='All IVs use market forward' if convention=='market' else 'Model IVs use own same-expiry forward'
            ax.set_title(f'VIX {days} days · {title}',loc='left',fontsize=11)
    handles,legend=axes[0,0].get_legend_handles_labels();fig.legend(handles,legend,loc='outside lower center',ncols=3,frameon=False)
    fig.suptitle('Forward convention changes the IV picture, not option prices\nSame strikes and prices in both columns; model-forward errors remain calibration targets.',fontsize=14)
    outputs=[]
    for suffix in ['png','svg']:
        path=directory/f'vix-forward-conventions.{suffix}';fig.savefig(path,dpi=140,bbox_inches='tight');outputs.append(path)
    plt.close(fig)
    manifest=dict(source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),directory/'iv_conventions.csv']},
        output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs},
        uncertainty='Only market bid/ask shading shown. No model-IV uncertainty band claimed: native IV also depends on estimated model forward and its covariance with price.')
    (directory/'plot_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--directory',type=Path,required=True)
    args=parser.parse_args();render(args.directory)
