import Mathlib.Tactic

/-! Finite algebraic bridges for a proposed exogenous stock-driver rotation.
These results do not prove a Brownian-process construction, VIX invariance,
continuous-time martingality, floating-point correctness, or a smile fit.
The recurrence theorem makes the numerical log-stock update an explicit premise;
source correspondence and joint Gaussian sampling need separate checks. -/
namespace M2Leverage

theorem rotation_coefficient_norm (theta r : ℝ)
    (ht0 : 0 ≤ theta) (ht1 : theta ≤ 1) (hr0 : 0 ≤ r) (hr1 : r ≤ 1) :
    (Real.sqrt (1-theta^2))^2 + (-theta*Real.sqrt (1-r))^2 +
      (-theta*Real.sqrt r)^2 = 1 := by
  have ht : 0 ≤ 1-theta^2 := by nlinarith
  simp only [mul_pow, neg_sq, Real.sq_sqrt ht,
    Real.sq_sqrt (sub_nonneg.mpr hr1), Real.sq_sqrt hr0]
  ring

theorem stress_covariance_contribution_nonpositive (theta r v ks kf : ℝ)
    (ht : 0 ≤ theta) (hs : 0 ≤ ks) (hf : 0 ≤ kf) :
    -theta * Real.sqrt v * (Real.sqrt (1-r)*ks+Real.sqrt r*kf) ≤ 0 := by
  have hsum : 0 ≤ Real.sqrt (1-r)*ks+Real.sqrt r*kf :=
    add_nonneg (mul_nonneg (Real.sqrt_nonneg _) hs)
      (mul_nonneg (Real.sqrt_nonneg _) hf)
  calc
    -theta * Real.sqrt v * (Real.sqrt (1-r)*ks+Real.sqrt r*kf) =
        -(theta * Real.sqrt v * (Real.sqrt (1-r)*ks+Real.sqrt r*kf)) := by ring
    _ ≤ 0 := neg_nonpos.mpr (mul_nonneg (mul_nonneg ht (Real.sqrt_nonneg _)) hsum)

theorem three_statistics_stock_reconstruction (theta : ℝ)
    (v dt dw db logStock : ℕ → ℝ) (h0 : logStock 0 = 0)
    (hstep : ∀ n, logStock (n+1) = logStock n - v n*dt n/2 +
      Real.sqrt (v n)*(Real.sqrt (1-theta^2)*dw n-theta*db n)) (n : ℕ) :
    logStock n = -(1/2 : ℝ)*(∑ i ∈ Finset.range n, v i*dt i) +
      Real.sqrt (1-theta^2)*(∑ i ∈ Finset.range n, Real.sqrt (v i)*dw i) -
      theta*(∑ i ∈ Finset.range n, Real.sqrt (v i)*db i) := by
  induction n with
  | zero => simp [h0]
  | succ n ih =>
    rw [hstep, ih]
    simp only [Finset.sum_range_succ]
    ring

#print axioms rotation_coefficient_norm
#print axioms stress_covariance_contribution_nonpositive
#print axioms three_statistics_stock_reconstruction
set_option pp.explicit true in
#check three_statistics_stock_reconstruction
end M2Leverage
