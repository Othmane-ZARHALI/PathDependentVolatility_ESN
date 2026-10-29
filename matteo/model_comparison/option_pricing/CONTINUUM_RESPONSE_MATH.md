# Continuous M2/stress response: implementation correspondence

For u in [0,1], define

    H_t(u) = G(X_t)^(1-u) / M(t,u)
             * exp(b*u*Z_t - .5*(b*u)^2*q(t)),
    M(t,u) = E_Q[G(X_t)^(1-u)].

The new OU starts at zero, is independent of the old Gaussian reservoir,
and has stationary variance one. q(t)=1-exp(-2*a*t).
The deterministic probability measure is Beta(m*kappa,(1-m)*kappa).
The floor/curve formula is

    v_t = f + (xi(t)-f) * integral H_t(u) mu(du).

These are research equations; there is no new Lean formalization for this family.

## Conditional pricing

At VIX expiry T and s>=T, independence gives

    E[H_s(u)|X_T,Z_T]
      = E[G(X_s)^(1-u)|X_T]/M(s,u)
        * exp(b*u*exp(-a*(s-T))*Z_T
               - .5*(b*u)^2*exp(-2*a*(s-T))*q(T)).

The variance curve multiplies this expression inside the future time integral.
The VIX square root is taken only after both integrals, then option payoffs.
No averaging of separately priced smiles occurs.

The implementation caches log(E[G^(1-u)|X_T]/M(s,u)) on a uniform u grid,
then uses cubic interpolation and positive Gauss-Jacobi quadrature for mu.
Cubic interpolation is an approximation and is tested off-grid; interpolating
logs preserves positivity, but does not itself enforce exact expected value one.
Normalizers are recomputed at each time node with a higher Gaussian rule.
All fitted xi levels remain free within their existing bounds.

## Positivity, normalization and interpretation

The old smooth cap bounds G above, and G is strictly positive. For fixed
u in [0,1], G^(1-u) <= max(1,G), so its normalizer exists, is finite and positive.
The Gaussian exponential factor has mean one. Independence, nonnegative
probability weights and Tonelli give E[integral H_t(u) mu(du)]=1.
Consequently E[v_t]=xi(t), and v_t>=f whenever xi(t)>=f.
These statements concern the mathematical model, not exact finite-quadrature
or finite-sample identities.

A fixed mu interval mass is the share of expected non-floor variance assigned
to that part of the response spectrum. It is not the probability of selecting
a regime. Realized component shares change with X and Z.

At u=0 the original normalized M2 response is recovered. At u=1 the response
is the normalized pure OU stress exponential. Intermediate responses retain
both factors, with weaker M2 exponent and stronger stress exposure as u grows.
As kappa approaches zero with m fixed, mu weakly approaches
(1-m)*delta_0 + m*delta_1. An ordinary Beta density with positive kappa has no
point mass at the endpoints; the boundary is reported separately.

## Forward constraint

Population forward: E[VIX_T]. Training equality: mean of the frozen weighted,
controlled VIX estimator equals its observed target for each of four expiries.
The optimizer enforces the latter to numerical tolerance. This is distinct from
E[VIX_T^2], whose population value is 10000/Delta times integral xi(s) ds.
Fresh sampling residuals and errors are needed before interpreting population
forward fit. No finite simulation alone certifies literal exact equality.
