"""
vectorized_sim.py
======================================
esn_base._sim_esn_with_params ALWAYS starts from r=0, z=0,
var_used_prev=0 -- it has no way to condition a simulation on an
arbitrary mid-path state. The VIX-smile test needs exactly that: an
"outer" ensemble of possible reservoir states at a future date N, and,
FROM EACH of those states, an "inner" ensemble of 30-day-forward paths
to price that day's own model-implied VIX.

This module reimplements the IDENTICAL equations (Eqs. 2-10 of the
companion note) in batched numpy, with 2 differences from
esn_base._sim_esn_with_params: (i) it accepts an arbitrary initial
state (r0, z0, var_used_prev0), and (ii) it simulates B independent
replicates at once instead of 1. Every equation is copied faithfully
from esn_base.py (same arch/matrices dict); the batching is a pure
vectorization, not a re-derivation. Validated against the original
scalar engine in validate_vectorized_sim.py.
"""
import numpy as np


def _softplus(x):
    return np.logaddexp(0.0, x)


def prepare(arch, mat, b0_delta=0.0):
    """Precompute everything that doesn't depend on the batch or the
    day loop (mirrors the per-call setup at the top of
    esn_base._sim_esn_with_params). b0_delta mirrors that function's
    own native b0_delta argument (additive offset to b0, Eq. 12)."""
    n_r, n_z = mat["n_r"], mat["n_z"]
    lam = mat["lambdas"]
    al = np.exp(-lam)                       # dt=1.0 throughout this note
    cl = np.sqrt(np.maximum(1.0 - al**2, 1e-14))
    az = mat["az"]
    azd = np.exp(-az)
    om = 1.0 - azd
    wz_vec = mat["z_readout_profile"] / np.sqrt(n_z)
    rc = arch["rough_orientation"] * mat["rough_scale"]
    m1 = mat["m1"] / n_r
    m2 = mat["m2"] / n_r
    q = mat["q"]
    fi, si = mat["fi"], mat["si"]
    zz, sz = mat["zz"], mat["sz"]
    jm_idx = np.clip(np.arange(n_z) - 1, 0, n_z - 1)
    jp_idx = np.clip(np.arange(n_z) + 1, 0, n_z - 1)
    return dict(n_r=n_r, n_z=n_z, al=al, cl=cl, azd=azd, om=om, wz_vec=wz_vec,
                rc=rc, m1=m1, m2=m2, q=q, fi=fi, si=si, zz=zz, sz=sz,
                jm_idx=jm_idx, jp_idx=jp_idx, b0=mat["b0"] + b0_delta,
                gamma_norm=arch["gamma_norm"], even_strength=arch["even_strength"],
                linear_strength=arch["linear_strength"], z_strength=arch["z_strength"],
                local_z_strength=arch["local_z_strength"], SIG_MIN=0.01 / np.sqrt(252.0))


def step(state, p, alpha, rng):
    """One vectorized day-step across a batch of B independent
    replicates. state = (r [B,Nr], z [B,Nz], var_used_prev [B])."""
    r, z, var_used_prev = state
    B = r.shape[0]
    eps = rng.standard_normal(B)

    qr = r @ p["q"]
    quad = p["m1"] * np.sum(r * r, axis=1) + p["m2"] * qr * qr
    eta = p["b0"] + p["rc"] * qr + z @ p["wz_vec"] + quad
    sp = _softplus(eta)
    var_raw = p["SIG_MIN"]**2 + sp**2
    var_used = (1 - alpha) * var_used_prev + alpha * var_raw
    sig_used = np.sqrt(np.maximum(var_used, 1e-30))

    dx = -0.5 * var_used + sig_used * eps          # dt=1, MU_S=0

    r_new = r * p["al"][None, :] + np.outer(eps, p["cl"])

    m_gate = np.maximum(0.0, 1.0 - p["gamma_norm"] * np.linalg.norm(r_new, axis=1) / np.sqrt(p["n_r"]))
    p1 = r_new[:, p["fi"]]
    p2 = r_new[:, p["si"]]
    ev = 0.7 * p1**2 + 0.3 * p2**2 - 1.0
    lin = 0.7 * p1 + 0.3 * p2
    lc = 0.5 * p["local_z_strength"] * (z[:, p["jm_idx"]] + z[:, p["jp_idx"]])
    u = p["sz"][None, :] * p1 * p2 + p["even_strength"] * ev - p["linear_strength"] * lin \
        + p["zz"][None, :] * z + lc
    z_new = p["azd"][None, :] * z + p["om"][None, :] * m_gate[:, None] * p["z_strength"] * np.tanh(u)

    return (r_new, z_new, var_used), dx, var_used


def simulate_batch(arch, mat, alpha, n_days, B, r0=None, z0=None, var_used_prev0=None, seed=0,
                    b0_delta=0.0):
    """B independent replicates, n_days forward, from a given (or
    default zero) initial state. Returns dx [n_days,B], var_used
    [n_days,B], and the final state tuple for further chaining."""
    p = prepare(arch, mat, b0_delta=b0_delta)
    r = np.zeros((B, p["n_r"])) if r0 is None else np.asarray(r0, float).copy()
    z = np.zeros((B, p["n_z"])) if z0 is None else np.asarray(z0, float).copy()
    vprev = np.zeros(B) if var_used_prev0 is None else np.asarray(var_used_prev0, float).copy()
    rng = np.random.default_rng(seed)

    dx_out = np.empty((n_days, B))
    v_out = np.empty((n_days, B))
    state = (r, z, vprev)
    for t in range(n_days):
        state, dx, v = step(state, p, alpha, rng)
        dx_out[t] = dx
        v_out[t] = v
    return dx_out, v_out, state
