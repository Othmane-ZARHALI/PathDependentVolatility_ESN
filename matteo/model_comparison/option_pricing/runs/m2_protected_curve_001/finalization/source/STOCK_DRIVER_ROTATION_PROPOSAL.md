# Proposed leverage experiment that preserves the variance-factor law

28 October 2026. **Proposal only; not implemented or calibrated.** Motivated by the protected nine-level SPX experiment. Its price loss improves on fresh paths, but the SPX downward skew remains missing. No model-family impossibility follows from that bounded search.

## A single risk-neutral leverage parameter

Let W 0 denote the old equity Brownian driver, which also drives the existing correlated M2 OU bank. Let Bs and Bf drive the slow/fast stress factors. The current model has the pair(Bs,Bf) independent of the entire old M2/ W 0 system. Preserve that joint factor law and all existing variance formulas.

Using the existing fixed fast-mixture weight r, define

$$dB^*=\sqrt{1-r}\,dB^s+\sqrt r\,dB^f,$$

$$dW^{stock}_\theta=\sqrt{1-\theta^2}\,dW^0-\theta\,dB^*,\qquad0\le\theta<1,$$

and price stocks through

$$dS_t/S_t=\sqrt{v_t}\,dW^{stock}_{\theta,t}.$$

Since the component Brownian drivers are independent and the squared coefficients sum to 1, Wstock is a standard Brownian driver. Positive theta makes positive stress shocks tend to accompany negative stock shocks. At theta=0 this is the current model. All maturities use the same parameter. The stress direction uses the existing r; it introduces no second fitted direction or factor.

This is a constructive way to maintain a valid joint covariance matrix. It differs from adding stress/equity correlations while demanding that all old equity/M2 correlations remain unchanged. Here

$$d\langle W^{stock}_\theta,B^s\rangle=-\theta\sqrt{1-r}\,dt,\qquad d\langle W^{stock}_\theta,B^f\rangle=-\theta\sqrt r\,dt,$$

and the old stock/M2 driver correlations are multiplied by sqrt(1-theta^2). The old-bank mutual correlations and all variance-factor dynamics stay fixed. Thus this is a change of the stock equation and its leverage interactions, not a claim that every previous stock correlation is preserved.

## Why VIX can remain protected

The entire joint law of the M2 states and both stress factors is unchanged. The variance v is the same function of those states and the same xi curve. None of their dynamics receives stock-return feedback. Therefore their conditional future variance law, and hence

$$VIX_T=100\sqrt{\frac1\Delta\int_T^{T+\Delta}E_Q[v_u\mid\mathcal F_T]du},$$

is unchanged. This preserves the VIX distribution and all its option prices/forwards, rather than merely constraining four mean values. It preserves existing numerical forward accuracy; it does not newly establish exact population-forward matching. This argument relies on exogenous factor dynamics and would need review for an endogenous return-feedback model.

This rotation therefore offers a more protective implementation of the previously discussed correlation idea. A different construction that alters dependence among M2/stress factors can change VIX and would need VIX refitting.

## What a correct implementation would require

An OU transition innovation is a weighted Brownian integral. Its normalized innovation is **not** the corresponding unweighted Brownian increment. For a stress rate a and step h, write

$$\epsilon=\frac{\sqrt{2a}\int_t^{t+h}e^{-a(t+h-u)}dB_u}{\sqrt{1-e^{-2ah}}}.$$

Its covariance with the Brownian increment is

$$c=Cov(\Delta B,\epsilon)=\frac{\sqrt{2a}(1-e^{-ah})}{a\sqrt{1-e^{-2ah}}}.$$

One can sample the correct joint pair using

$$\Delta B=c\epsilon+\sqrt{h-c^2}\,\eta,$$

with a new independent standard normal eta. Small-step arithmetic needs a stable evaluation and a covariance audit. Then form the rotated stock increment from the old W 0 increment and both correct stress Brownian increments. Reusing Gaussian state draws is possible, but stock paths must be rebuilt. Simply mixing normalized OU innovations into the equity increments would change the intended joint transition law.

The present independent-stress Black-delta controls must **not** be reused unchanged: after rotation, conditioning on future stress also conditions on part of the stock driver. The zero-mean hedge and conditional-forward arguments change. A new reference/conditional pricing derivation and independent covariance and price checks are needed before fitting theta.

This proposal fixes the risk-neutral factor law. The physical-measure interpretation of the old parameter named equity risk premium would need explicit review after changing the stock driver; no physical-recalibration claim is made here. Stock-martingale and numerical accuracy checks also remain separate from the Brownian covariance construction.
