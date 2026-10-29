# Numerical forward repair, 27 October 2026

The initial continuous experiment002 met four training-estimator equalities but
fresh55-day residual was -.1656 versus .0248SE. The cheap two-atom reference
refitted on the same8192points failed similarly; this is not evidence of a
continuous-model error. Training forward SEs were .038–.098 and old frozen
controls were noisy (55day training control mean -.195). Reject that calibration
as a reliable fulfillment of the forward condition; preserve all results.

Keep the mathematical model, Beta family, bounds, all217native quotes, four
forwards, fitted floor and seven independent active levels unchanged.
Use the exact expectation identity
 E[P_cont]=E[P_boundary]+E[P_cont-P_boundary]
for forwards and option payoffs. The boundary is the matching kappa=0 model
with the same b,m,half-life,floor and curve. Estimate its first term on32768
existing parent checking paths (now training); estimate the paired difference
on the8192 powered training paths already computed. These are nested nets;
eight-net SEs combine whole matching net blocks, retaining covariance.
The identity changes only numerical variance reduction, not the model.

Three starts: prior endpoint winner, prior best kappa>=.1 solution, and prior
secondary local minimum. SLSQP four equality constraints, same eligibility1e-5.
Select main by training option loss and secondary by training loss amongkappa>=.1.
No fitting or selection on checking results.

Checking uses new65536points/expiry for the cheap boundary term (8x8192 new
scrambles, seed20261027117), plus the existing16384 continuous-minus-boundary
checking paths. The difference stream is explicitly reused after investigation;
this is a mixed numerical check, not a wholly independent fresh study. Report
sampling SEs, training uncertainty, paired differences and variance diagnostics.
Original experimental fit002 and original production baseline remain unchanged.
No SPX pricing, new model factor, curve smoothing, clipping or new Lean theorem.
