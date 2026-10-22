"""Fixed component allocations for a frozen defensive Gaussian mixture."""
from dataclasses import dataclass
import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc


@dataclass(frozen=True)
class StratifiedDraw:
    normals: np.ndarray
    likelihood_weights: np.ndarray
    allocation_factors: np.ndarray
    components: np.ndarray
    component_points: np.ndarray
    power: int
    owner: object

    @property
    def weights(self):
        """Weights for a plain mean; includes the exact allocation correction."""
        return self.likelihood_weights*self.allocation_factors


class StratifiedGaussianMixture:
    """Estimate sum_k pi_k mean_k[(f/g) H] with fixed positive counts.

    The proposal g=sum pi_k phi_k is unchanged. For N total draws and n_k from
    component k, the plain-mean coefficient is (N*pi_k/n_k)*(f/g). This is
    deterministic allocation correction, never sample-weight normalization.
    Complete independent sets of component scrambles are uncertainty units.
    """
    def __init__(self, mixture, base_counts=(32,6,6,5,5,5,5), base_power=6):
        counts=np.asarray(base_counts)
        if (not isinstance(base_power,int) or base_power<1 or counts.shape!=mixture.probabilities.shape
                or not np.issubdtype(counts.dtype,np.integer) or np.any(counts<=0)
                or int(counts.sum())!=2**base_power):
            raise ValueError('Require positive integer component counts summing to 2**base_power')
        self.mixture,self.base_power=mixture,base_power
        self.base_counts=counts.astype(int,copy=True)
        self.base_counts.setflags(write=False)

    def counts(self,power):
        if not isinstance(power,int) or power<self.base_power:
            raise ValueError('Power below declared allocation base')
        return self.base_counts*2**(power-self.base_power)

    def draw(self,power,seed):
        counts=self.counts(power)
        normals,components,points,factors=[],[],[],[]
        for k,n in enumerate(counts):
            component_seed=int(np.random.SeedSequence([seed,k]).generate_state(1)[0])
            # A non-power-of-two allocation uses a prefix of the next full net.
            # No assertion of full Sobol balance is made for that prefix.
            u=qmc.Sobol(self.mixture.shifts.shape[1],scramble=True,bits=52,rng=component_seed).random_base2(int(n-1).bit_length())[:n]
            if np.any((u<=0)|(u>=1)):
                raise FloatingPointError('Sobol endpoint')
            normals.append(ndtri(u)+self.mixture.shifts[k])
            components.append(np.full(n,k,int));points.append(np.arange(n))
            factors.append(np.full(n,2**power*self.mixture.probabilities[k]/n))
        z=np.concatenate(normals)
        arrays=[z,self.mixture.weights(z),np.concatenate(factors),np.concatenate(components),np.concatenate(points)]
        for a in arrays:
            a.setflags(write=False)
        return StratifiedDraw(*arrays,power,self)

    def prefix_mask(self,draw,power):
        if draw.owner is not self or power>draw.power:
            raise ValueError('Prefix belongs to a different sampler or exceeds fine draw')
        counts=self.counts(power)
        return draw.component_points<counts[draw.components]
