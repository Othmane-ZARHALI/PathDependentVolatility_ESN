"""Single-call convention and deterministic, continuous nested path budgets.

Prices are undiscounted. Market conversion uses the market forward; model put
reconstruction must use the model forward, particularly for VIX.
"""
import numpy as np
from scipy.special import ndtr

TIERS = (1024, 2048, 4096)


def canonical_market(quotes):
    q = quotes.copy()
    adjustment = np.where(q.option_type.to_numpy() == 'put',
                          q.forward.to_numpy()-q.strike.to_numpy(), 0.)
    for column in ['price_bid_U', 'price_mid_U', 'price_ask_U']:
        q[column] = q[column].to_numpy()+adjustment
    q['original_option_type'] = q.option_type
    q['option_type'] = 'call'
    return q


def tier_weights(distance):
    d = np.asarray(distance, float)
    if not np.isfinite(d).all() or np.any(d < 0):
        raise ValueError('Require finite nonnegative standardized distances')
    a = np.clip(d-1, 0, 1); a = a*a*(3-2*a)
    b = np.clip(d-2, 0, 1); b = b*b*(3-2*b)
    return np.stack([1-a, a-b, b], axis=-1)


def sample_coefficients(distance, tiers=TIERS):
    if len(tiers) != 3 or any(int(n) != n or n < 2 for n in tiers) or not np.all(np.diff(tiers) > 0):
        raise ValueError('Require three increasing integer budgets')
    weights = tier_weights(distance).reshape(-1, 3)
    a = np.zeros((tiers[-1], len(weights)))
    for j, n in enumerate(tiers):
        a[:n] += weights[:, j]/n
    return a


def iid_statistics(draws, distance, tiers=TIERS):
    """Unbiased variance estimate for a deterministic mixture of iid prefixes.

    Each column uses identical-law independent draws, despite shared prefixes.
    Var(sum a_i Z_i)=Var(Z)*sum a_i^2. Sample variance uses the active prefix.
    Cross-strike covariance is not represented by the pointwise SEs.
    """
    y = np.asarray(draws, float)
    a = sample_coefficients(distance, tiers)
    if y.shape != a.shape or not np.isfinite(y).all():
        raise ValueError('Draws must have maximum-budget by quote shape')
    counts = np.sum(a > 0, axis=0)
    variance = np.array([np.var(y[:n, j], ddof=1) for j, n in enumerate(counts)])
    squares = np.sum(a*a, axis=0)
    return dict(price=np.sum(a*y, axis=0), se=np.sqrt(variance*squares),
                actual_paths=counts, effective_paths=1/squares, coefficients=a)


def scrambled_statistics(draws, distance, tiers=TIERS):
    """Each scramble is one replication, including all shared-prefix blending."""
    y = np.asarray(draws, float)
    if y.ndim != 3 or len(y) < 2 or any(n % len(y) for n in tiers):
        raise ValueError('Require independent nets and budgets divisible by net count')
    sizes = tuple(n//len(y) for n in tiers)
    a = sample_coefficients(distance, sizes)
    if y.shape[1:] != a.shape:
        raise ValueError('Net/point/quote shape inconsistent with budgets')
    b = np.einsum('nq,rnq->rq', a, y)
    return dict(price=b.mean(0), se=b.std(0, ddof=1)/np.sqrt(len(b)),
                actual_paths=len(b)*np.sum(a > 0, axis=0),
                effective_paths=len(b)/np.sum(a*a, axis=0), blocks=b,
                coefficients=np.tile(a/len(b), (len(b), 1)))


def reference_calls(r, dr, ref, dref, weights, strikes, total_variance, dvariance):
    """Known lognormal call plus zero-mean differences; analytic Jacobian.

    All strikes have forward 1 here. No put branch or strike-anchor switching.
    dvariance has one entry per model parameter. Ref is driven by the same
    equity increments with deterministic variance xi(t); it is NOT normalized
    by a sample moment. The derivative includes the reference and control.
    """
    k = np.asarray(strikes, float); a = float(total_variance)
    if a <= 0 or np.any(k <= 0):
        raise ValueError('Positive reference variance and strikes required')
    d1 = -np.log(k)/np.sqrt(a)+.5*np.sqrt(a)
    delta = ndtr(d1); phi = np.exp(-.5*d1*d1)/np.sqrt(2*np.pi)
    known = delta-k*ndtr(d1-np.sqrt(a))
    dknown = phi/(2*np.sqrt(a))
    ddelta = phi*(np.log(k)/(2*a**1.5)+.25/np.sqrt(a))
    actual = np.maximum(r[:, None]-k, 0.)
    reference = np.maximum(ref[:, None]-k, 0.)
    draws = known+weights[:, None]*(actual-reference-delta*(r-ref)[:, None])
    derivative = (dknown[:, None]*dvariance[None, :])[None, :, :]+weights[:, None, None]*(
        ((r[:, None] > k)-delta)[..., None]*dr[:, None, :]
        -((ref[:, None] > k)-delta)[..., None]*dref[:, None, :]
        -(r-ref)[:, None, None]*ddelta[None, :, None]*dvariance[None, None, :])
    return draws, derivative


def gaussian_controls(z, weights, dimensions=6):
    z = np.asarray(z)[:, :dimensions]; w = np.asarray(weights)
    # z has target N(0,I), even though proposal draws are shifted.
    return np.column_stack([w-1, w[:, None]*z, w[:, None]*(z*z-1)])


def fit_control_coefficients(control, payoffs):
    c, y = np.asarray(control), np.asarray(payoffs)
    # Independently generated pilot only. Centering is for OLS coefficients;
    # evaluation uses the theoretical zero means, never pilot/sample means.
    c = c-c.mean(0); y = y-y.mean(0)
    return np.linalg.lstsq(c, y, rcond=1e-10)[0]
