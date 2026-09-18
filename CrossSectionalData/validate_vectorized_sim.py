import numpy as np
import esn_base as E
from joint_calibration import OPT_ARCH, P_FIXED
from vectorized_sim import simulate_batch

SP500_PARAMS = dict(H=0.386, rough_scale=0.672, lam_lo=0.00001,
                     lam_hi=4.8809, m1=0.1745, alpha=0.220)
ip = dict(P_FIXED, H=SP500_PARAMS['H'], rough_scale=SP500_PARAMS['rough_scale'],
          lam_lo=SP500_PARAMS['lam_lo'], lam_hi=SP500_PARAMS['lam_hi'], m1=SP500_PARAMS['m1'])
mat = E.build_esn_matrices(OPT_ARCH, ip)
alpha = SP500_PARAMS['alpha']

T = 3000
N_scalar = 400

# ---- scalar engine, N_scalar independent runs (different seeds) ----
xs_scalar = []
vs_scalar = []
for k in range(N_scalar):
    x, v, vr = E._sim_esn_with_params(9_000_000 + k, T, OPT_ARCH, mat, alpha=alpha)
    xs_scalar.append(x)
    vs_scalar.append(v)
xs_scalar = np.array(xs_scalar)   # (N_scalar, T)
vs_scalar = np.array(vs_scalar)

# ---- vectorized batch engine, same T, N_scalar replicates, one call ----
dx_vec, v_vec, _ = simulate_batch(OPT_ARCH, mat, alpha, T, N_scalar, seed=9_500_000)
dx_vec = dx_vec.T   # (N_scalar, T)
v_vec = v_vec.T

print("=== distributional agreement check (independent seeds both sides) ===")
print(f"{'stat':<20}{'scalar engine':>16}{'vectorized':>16}")
print(f"{'mean(x)':<20}{xs_scalar.mean():>16.6f}{dx_vec.mean():>16.6f}")
print(f"{'std(x)':<20}{xs_scalar.std():>16.6f}{dx_vec.std():>16.6f}")
print(f"{'skew(x)':<20}{((xs_scalar-xs_scalar.mean())**3).mean()/xs_scalar.std()**3:>16.4f}"
      f"{((dx_vec-dx_vec.mean())**3).mean()/dx_vec.std()**3:>16.4f}")
c = xs_scalar - xs_scalar.mean(); kv = (c**4).mean()/xs_scalar.var()**2 - 3
c2 = dx_vec - dx_vec.mean(); kv2 = (c2**4).mean()/dx_vec.var()**2 - 3
print(f"{'excess kurtosis':<20}{kv:>16.4f}{kv2:>16.4f}")
print(f"{'mean(v)':<20}{vs_scalar.mean():>16.8f}{v_vec.mean():>16.8f}")
print(f"{'std(v)':<20}{vs_scalar.std():>16.8f}{v_vec.std():>16.8f}")

# acf(|x|) lag 1,5,10 as a dynamics check (not just marginal moments)
def acf_pooled(X, lag):
    vals = []
    for row in X:
        a = row - row.mean()
        vals.append(np.mean(a[:-lag]*a[lag:])/np.var(a))
    return np.mean(vals)

for lag in [1,5,10]:
    a1 = acf_pooled(np.abs(xs_scalar), lag)
    a2 = acf_pooled(np.abs(dx_vec), lag)
    print(f"acf(|x|) lag {lag:2d}:  scalar={a1:.4f}   vectorized={a2:.4f}")
