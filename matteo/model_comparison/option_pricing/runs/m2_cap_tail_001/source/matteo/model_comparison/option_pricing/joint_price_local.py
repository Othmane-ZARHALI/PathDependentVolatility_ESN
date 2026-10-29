"""Cheap affine SPX AND VIX price/forward model at a true simulated anchor.

No simulator/engine is held here. The original interval loss is evaluated on
approximate prices; actual nonlinear prices alone can accept a proposal.
Model VIX forward derivatives are included in native put-call parity.
"""
import numpy as np
from .all_maturity_core import InSampleLoss


class JointPriceLocal:
    def __init__(self, parameters, anchor, canonical_quotes, native_quotes, forwards):
        self.parameters = np.asarray(parameters, float).copy()
        self.anchor = anchor
        self.loss = InSampleLoss(canonical_quotes, forwards)
        self.native_loss = InSampleLoss(native_quotes, forwards)
        self.quotes = self.loss.quotes
        self.put = native_quotes.option_type.eq('put').to_numpy()
        self.forward_index = np.full(len(native_quotes), -1, dtype=int)
        for j, day in enumerate(forwards.expiry_days):
            mask = native_quotes.instrument.eq('VIX') & native_quotes.expiry_days.eq(day)
            self.forward_index[mask.to_numpy()] = j
        assert np.all(self.forward_index[native_quotes.instrument.eq('VIX')] >= 0)
        self.calls = 0

    def evaluate(self, parameters):
        self.calls += 1
        p = np.asarray(parameters, float)
        delta = p - self.parameters
        pj = self.anchor['price_jacobian']
        fj = self.anchor['forward_jacobian']
        prices = self.anchor['prices'] + pj @ delta
        forwards = self.anchor['forwards'] + fj @ delta
        model_forward = self.quotes.forward.to_numpy().copy()
        forward_derivative = np.zeros_like(pj)
        vix = self.forward_index >= 0
        model_forward[vix] = forwards[self.forward_index[vix]]
        forward_derivative[vix] = fj[self.forward_index[vix]]
        native = prices - self.put * (model_forward - self.quotes.strike.to_numpy())
        nj = pj - self.put[:, None] * forward_derivative
        results = []
        gradients = []
        vix_gradients = []
        penalty_jac = np.zeros((2, 11))
        penalty_jac[:, 7:10] = np.sqrt(.1) * np.array([[-1, 1, 0], [0, -1, 1]])
        for loss, value, derivative in [(self.loss, prices, pj), (self.native_loss, native, nj)]:
            r = loss.evaluate(value, forwards, p[7:10])
            q = loss.quotes
            outside = (value < q.price_bid_U.to_numpy()) | (value > q.price_ask_U.to_numpy())
            residual_jac = loss.sqrtweights[:, None] * outside[:, None] * derivative / q.price_scale_U.to_numpy()[:, None]
            full_jac = np.concatenate([residual_jac, fj / (.5 * np.sqrt(len(forwards))), penalty_jac])
            gradients.append(2 * r['residuals'] @ full_jac)
            is_vix = q.instrument.eq('VIX').to_numpy()
            vr = r['residuals'][:len(q)][is_vix]
            fr = r['residuals'][len(q):len(q)+len(forwards)]
            vix_gradients.append(2 * vr @ residual_jac[is_vix] + 2 * fr @ (fj / (.5 * np.sqrt(len(forwards)))))
            r['total'] += .1 / 3
            r['components']['state_penalty'] = .1 / 3
            results.append(r)
        common, native_result = results
        cv = common['components']['VIX'] + common['components']['forwards']
        nv = native_result['components']['VIX'] + native_result['components']['forwards']
        return dict(criteria=np.array([common['total'], native_result['total'], cv, nv]),
                    criteria_gradients=np.stack([*gradients, *vix_gradients]),
                    merit=nv, merit_gradient=vix_gradients[1],
                    prices=prices, native_prices=native, forwards=forwards,
                    canonical=common, native=native_result,
                    approximate_prices=True, simulator_calls=0)
