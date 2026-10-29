import Mathlib.MeasureTheory.Integral.Bochner.Basic
import Mathlib.Tactic

/-! Restricted mechanism checks, not a theorem about M2's implied smile.
A,B represent positive exponential branch values; a,b their log sensitivities.
Identifying the expression with a log derivative is an analytic correspondence
obligation. The portfolio theorem assumes the entire payoff inequality. -/
open MeasureTheory
namespace M2Wing

theorem response_curvature_numerator (A B a b : ℝ) :
    (a^2*A+b^2*B)*(A+B)-(a*A+b*B)^2 = A*B*(a-b)^2 := by
  ring

theorem response_curvature_positive (A B a b : ℝ)
    (hA : 0 < A) (hB : 0 < B) (hab : a ≠ b) :
    0 < ((a^2*A+b^2*B)*(A+B)-(a*A+b*B)^2)/(A+B)^2 := by
  rw [response_curvature_numerator]
  exact div_pos (mul_pos (mul_pos hA hB) (sq_pos_of_ne_zero (sub_ne_zero.mpr hab)))
    (sq_pos_of_pos (add_pos hA hB))

theorem call_second_moment_payoff_bound (y K : ℝ) (hy : 0 ≤ y) (hK : 0 ≤ K) :
    4*K*max (y-K) 0 ≤ y^2 := by
  by_cases h : K ≤ y
  · rw [max_eq_left (sub_nonneg.mpr h)]
    nlinarith [sq_nonneg (y-2*K)]
  · rw [max_eq_right (by linarith : y-K ≤ 0)]
    nlinarith

theorem portfolio_moment_obstruction {Ω : Type*} [MeasurableSpace Ω]
    (μ : Measure Ω) (X P : Ω → ℝ) (budget cost : ℝ)
    (hX : Integrable (fun ω => (X ω)^2) μ) (hP : Integrable P μ)
    (hdom : ∀ ω, P ω ≤ (X ω)^2)
    (hcost : cost ≤ ∫ ω, P ω ∂μ)
    (hbudget : (∫ ω, (X ω)^2 ∂μ) = budget) : cost ≤ budget := by
  rw [← hbudget]
  exact le_trans hcost (integral_mono hP hX hdom)

#print axioms response_curvature_numerator
#print axioms response_curvature_positive
#print axioms call_second_moment_payoff_bound
#print axioms portfolio_moment_obstruction
set_option pp.explicit true in
#check portfolio_moment_obstruction
end M2Wing
