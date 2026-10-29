"""Independent-level MLMC statistics and bounded integer variance/cost allocation."""
import numpy as np


def allocate_levels(variances, costs, budget, minimum=64, quantum=64):
    """Minimize sum(V_l/N_l) with measured cost sum(C_l*N_l)<=budget.

    Continuous allocation is proportional to sqrt(V/C), with active minimum
    constraints. Round down to path blocks, then spend affordable blocks by
    marginal variance reduction per cost. Inputs must come from design data,
    independent of final checking draws. Empirical zero variance is not an
    accuracy certificate; the minimum count still applies to every level.
    """
    v,c=np.asarray(variances,float),np.asarray(costs,float)
    if (v.ndim!=1 or not len(v) or c.shape!=v.shape or not np.isfinite(v).all() or
            not np.isfinite(c).all() or np.any(v<0) or np.any(c<=0) or
            not np.isfinite(budget) or budget<=0 or not isinstance(minimum,int) or
            not isinstance(quantum,int) or quantum<1 or minimum<quantum or minimum%quantum):
        raise ValueError('Invalid allocation inputs')
    if minimum*c.sum()>budget:
        raise ValueError('Budget cannot fund every minimum level count')
    n=np.full(len(v),float(minimum));free=v>0
    while np.any(free):
        remaining=budget-float(n[~free]@c[~free])
        candidate=remaining*np.sqrt(v[free]/c[free])/np.sum(np.sqrt(v[free]*c[free]))
        indices=np.flatnonzero(free);below=candidate<minimum
        if not below.any():
            n[indices]=candidate;break
        free[indices[below]]=False
    n=np.maximum(minimum,np.floor(n/quantum).astype(int)*quantum)
    while True:
        affordable=(n@c+quantum*c)<=budget+1e-12*budget
        gain=np.where(affordable,v*(1/n-1/(n+quantum))/(quantum*c),-1.)
        best=int(np.argmax(gain))
        if gain[best]<=0:
            break
        n[best]+=quantum
    return n


def independent_level_statistics(level_draws):
    """Sum independent level means; SE²=sum unbiased sample variances/N_l.

    Each array's first axis is iid paths within that level. Level samples must
    be mutually independent conditional on frozen pilots/allocation. Supports
    vector payoffs; this does not bound finest-grid bias or unresolved tails.
    """
    if not level_draws:
        raise ValueError('Need at least one level')
    arrays=[np.asarray(x,float) for x in level_draws]
    if any(x.ndim<1 or len(x)<2 or x.shape[1:]!=arrays[0].shape[1:] or not np.isfinite(x).all() for x in arrays):
        raise ValueError('Require compatible finite arrays with at least two paths')
    return (sum(x.mean(axis=0) for x in arrays),
            np.sqrt(sum(x.var(axis=0,ddof=1)/len(x) for x in arrays)))
