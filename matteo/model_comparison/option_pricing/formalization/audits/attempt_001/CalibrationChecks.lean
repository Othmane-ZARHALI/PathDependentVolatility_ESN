import Mathlib.Probability.Distributions.Gaussian.Real
import Mathlib.MeasureTheory.Integral.Bochner.Basic
import Mathlib.Tactic

/-!
Scoped checks for fixed-model curve calibration. These results do not certify
floating-point Python, Monte Carlo tails, fitted parameters or a market fit.
The Gaussian one-step law is proved using Mathlib's Gaussian MGF. The control
reuse result assumes integrable, known-zero controls under ONE FIXED measure.
Those are explicit correspondence obligations for the numerical implementation.
-/

open MeasureTheory ProbabilityTheory
open scoped BigOperators

namespace M2Calibration

theorem gaussian_exponential_mean_one (v : ℝ≥0) :
    (∫ z : ℝ, Real.exp (z - (v : ℝ) / 2) ∂gaussianReal 0 v) = 1 := by
  have hmgf : (∫ z : ℝ, Real.exp z ∂gaussianReal 0 v) = Real.exp ((v : ℝ) / 2) := by
    simpa [mgf] using congrFun (mgf_fun_id_gaussianReal (μ := 0) (v := v)) (1 : ℝ)
  simp_rw [Real.exp_sub]
  rw [integral_div, hmgf, div_self (ne_of_gt (Real.exp_pos _))]

theorem gaussian_stock_step_mean (v : ℝ≥0) (stock : ℝ) :
    (∫ z : ℝ, stock * Real.exp (z - (v : ℝ) / 2) ∂gaussianReal 0 v) = stock := by
  rw [integral_const_mul, gaussian_exponential_mean_one, mul_one]

theorem gaussian_hedge_step_mean_zero (v : ℝ≥0) (stock delta : ℝ) :
    (∫ z : ℝ, delta * (stock * Real.exp (z - (v : ℝ) / 2) - stock)
      ∂gaussianReal 0 v) = 0 := by
  have hi : Integrable (fun z : ℝ => Real.exp z) (gaussianReal 0 v) := by
    simpa using integrable_exp_mul_gaussianReal (μ := 0) (v := v) (1 : ℝ)
  have hshift : Integrable (fun z : ℝ => Real.exp (z - (v : ℝ) / 2))
      (gaussianReal 0 v) := by
    simp_rw [Real.exp_sub]
    exact hi.div_const _
  rw [integral_const_mul, integral_sub (hshift.const_mul stock) (integrable_const stock),
    gaussian_stock_step_mean]
  simp

theorem frozen_anchor_control_reuse {Ω Θ : Type*} [MeasurableSpace Ω]
    (μ : Measure Ω) (n : ℕ) (H : Θ → Ω → ℝ) (C : Fin n → Ω → ℝ)
    (beta : Fin n → ℝ) (hH : ∀ θ, Integrable (H θ) μ)
    (hC : ∀ i, Integrable (C i) μ) (hzero : ∀ i, ∫ ω, C i ω ∂μ = 0) (θ : Θ) :
    (∫ ω, (H θ ω - ∑ i, beta i * C i ω) ∂μ) = ∫ ω, H θ ω ∂μ := by
  have hint : Integrable (fun ω => ∑ i, beta i * C i ω) μ := by
    exact integrable_finsetSum Finset.univ (fun i _ => (hC i).const_mul (beta i))
  rw [integral_sub (hH θ) hint, integral_finsetSum]
  · simp_rw [integral_const_mul, hzero, mul_zero]
    simp
  · intro i _
    exact (hC i).const_mul (beta i)

theorem finite_importance_identity (n : ℕ) (p q H : Fin n → ℝ)
    (hq : ∀ i, q i ≠ 0) :
    (∑ i, q i * (p i / q i) * H i) = ∑ i, p i * H i := by
  apply Finset.sum_congr rfl
  intro i _
  field_simp

def integratedVariance (n : ℕ) (floorPart floor : ℝ) (A xi : Fin n → ℝ) : ℝ :=
  floorPart + ∑ i, A i * (xi i - floor)

theorem integrated_variance_affine (n : ℕ) (floorPart floor t : ℝ)
    (A xi eta : Fin n → ℝ) :
    integratedVariance n floorPart floor A (fun i => (1-t)*xi i+t*eta i) =
      (1-t)*integratedVariance n floorPart floor A xi +
      t*integratedVariance n floorPart floor A eta := by
  unfold integratedVariance
  have hterm : ∀ i, A i * ((1-t)*xi i+t*eta i-floor) =
      (1-t)*(A i*(xi i-floor))+t*(A i*(eta i-floor)) := by
    intro i
    ring
  simp_rw [hterm, Finset.sum_add_distrib, ← Finset.mul_sum]
  ring

theorem integrated_variance_positive (n : ℕ) (floorPart floor : ℝ)
    (A xi : Fin n → ℝ) (hc : 0 < floorPart) (ha : ∀ i, 0 ≤ A i)
    (hx : ∀ i, floor ≤ xi i) :
    0 < integratedVariance n floorPart floor A xi := by
  unfold integratedVariance
  have hs : 0 ≤ ∑ i, A i * (xi i-floor) := by
    exact Finset.sum_nonneg (fun i _ => mul_nonneg (ha i) (sub_nonneg.mpr (hx i)))
  linarith

theorem weighted_residual_norm_identity (n : ℕ) (w r : Fin n → ℝ)
    (hw : ∀ i, 0 ≤ w i) :
    (∑ i, (Real.sqrt (w i) * r i)^2) = ∑ i, w i * (r i)^2 := by
  apply Finset.sum_congr rfl
  intro i _
  rw [mul_pow, Real.sq_sqrt (hw i)]

theorem penalty_residual_identity (a b c : ℝ) :
    (Real.sqrt (1/10 : ℝ) * (b-a))^2 + (Real.sqrt (1/10 : ℝ) * (c-b))^2 =
      (1/10 : ℝ) * ((b-a)^2+(c-b)^2) := by
  rw [mul_pow, mul_pow, Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 1/10)]
  ring

theorem finite_noise_square_decomposition (n : ℕ) (p noise : Fin n → ℝ) (a : ℝ)
    (hmass : ∑ i, p i = 1) (hmean : ∑ i, p i * noise i = 0) :
    (∑ i, p i * (a+noise i)^2) = a^2 + ∑ i, p i * (noise i)^2 := by
  have hterm : ∀ i, p i*(a+noise i)^2 = p i*a^2+2*a*(p i*noise i)+p i*(noise i)^2 := by
    intro i
    ring
  simp_rw [hterm, Finset.sum_add_distrib, Finset.sum_mul, ← Finset.mul_sum]
  rw [hmass, hmean]
  ring

#print axioms gaussian_exponential_mean_one
#print axioms gaussian_stock_step_mean
#print axioms gaussian_hedge_step_mean_zero
#print axioms frozen_anchor_control_reuse
#print axioms finite_importance_identity
#print axioms integrated_variance_affine
#print axioms integrated_variance_positive
#print axioms weighted_residual_norm_identity
#print axioms penalty_residual_identity
#print axioms finite_noise_square_decomposition

set_option pp.explicit true in
#check frozen_anchor_control_reuse
set_option pp.explicit true in
#check gaussian_hedge_step_mean_zero
set_option pp.explicit true in
#check integrated_variance_affine
set_option pp.explicit true in
#check finite_noise_square_decomposition

end M2Calibration
