"""Native option IV comparison using each model's own checking forward."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .plot_joint_curve_fit import implied_vol


def render(run,out):
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    for ax,day in zip(axes.flat,[27,55,118,153]):
        q=pd.read_csv(run/'fresh_check/native_iv_bank.csv').query('expiry_days == @day')
        def market_iv(col):return np.array([implied_vol(getattr(a,col),a.forward,a.strike,a.T_years,a.option_type) for a in q.itertuples()],float)
        ax.fill_between(q.K_over_F,market_iv('price_bid_U'),market_iv('price_ask_U'),color='.3',alpha=.18,label='Market bid–ask')
        ax.plot(q.K_over_F,q.market_native_IV,'k-',label='Market mid')
        for name,color,ls,label in [('two_speed','#236eb0','--','Previous fitted two-speed model'),('bank','#d04a35','-','Fixed geometric bank')]:
            z=pd.read_csv(run/'fresh_check'/f'native_iv_{name}.csv').query('expiry_days == @day');y=z.model_native_IV
            ax.plot(z.K_over_F,y,color=color,ls=ls,lw=1.8,label=label);invalid=~np.isfinite(y)
            if invalid.any():ax.plot(z.K_over_F[invalid],np.full(invalid.sum(),.02),'|',color=color,transform=ax.get_xaxis_transform())
        ax.set_xscale('log');ax.set(title=f'VIX {day} days',xlabel='Strike / market forward',ylabel='Native Black implied volatility (%)');ax.grid(alpha=.18)
    handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncols=2,frameon=False)
    fig.suptitle('Native puts and calls · each model uses its own checking forward\nFour fixed half-lives: 7, 28, 112, 448 days · unavailable IVs marked at foot',fontsize=12)
    out.mkdir(parents=True,exist_ok=True)
    for ext in ['png','svg']:fig.savefig(out/f'vix-stress-bank-native-IV.{ext}',dpi=155,bbox_inches='tight')
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();render(a.run,a.output)
