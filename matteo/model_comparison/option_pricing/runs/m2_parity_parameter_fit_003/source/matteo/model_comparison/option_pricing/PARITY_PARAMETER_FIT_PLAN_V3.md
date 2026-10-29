# Full validation of the near-feasible VIX proposal

Run002 completed its first20-evaluation constrained local search but had no point
with reserve slack>=-1e-8. Its last point has minimum reserve slack -2.93e-5:
the actual economic loss criteria all improve, and the optional reserve margin
is missed by only about0.003 loss units. This is far smaller than the reserve
itself. Do not spend another full solve merely to close that numerical margin.

Run003 reuses the EXACT saved/validated full4096 baseline prices and Jacobians
from002; sources, checkpoint files and the chosen warm local proposal are frozen.
No model/pricer/true-Jacobian change. The warm point is selected solely from
training local evaluations with reserve slack>=-1e-4. It receives immediate FULL
true4096 repricing/Jacobian. Only the original four strict true-price acceptance
criteria can accept it. The numerical reserve is guidance, not the final test.

Then at most2 further VIX-focused constrained SLSQP cycles, same settings as V2,
with reserve feasibility tolerance1e-4 normalized (versus2%/0.5% reserves), and
unchanged strict true-economic guards. Maximum75 optimization minutes including
the warm full repricing. All11 parameters/bounds/beta_J>0, geometry/Q/state,
OTM input scope4096 samples/GH/time settings remain unchanged.

Full checkpoint replay avoids another identical SPX simulation. Reuse the
production Jacobian screens already passed in001/002. Save all new true anchor
Jacobians. Fresh checking input hashes are already frozen; checking prices
remain unread until final selection. Preserve prototype001, near-feasible002,
and rejected/accepted true proposals without relabelling them as completed fits.

After the warm proposal, prioritize the ORIGINAL OTM VIX-plus-forward loss in
the guidance objective, while still constraining both converted/native VIX and
joint criteria. This aligns the search with actual quoted OTM contracts; price
simulation still uses one call payoff. The native-priority merit and gradient
are exactly the fourth criterion already covered by the production/all-column
constraint derivative tests. Reported objectives/weights and all four strict
true-price acceptance conditions are unchanged. This avoids spending the VIX
improvement mostly on the synthetic convention while native improvement sits
at its smallest allowed margin.
