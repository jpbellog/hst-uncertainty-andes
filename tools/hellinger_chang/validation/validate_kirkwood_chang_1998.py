"""
Validation of combine.welch_df / welch_critical_value / conservative_critical_value
against Kirkwood and Chang (1998).

Targets from the paper:
- Full data (Appendix B, p. 99): the example gives multiplier d = 3.049 and v = 42.34.
  The form is v = (12*Var + 2*E^2)/(3*Var - 2*E^2), d = E*(1 - 2/v), with E = 3+2a,
  Var = 6+14a+2c (b=0 for rotations, Theorem 5.1.2). Reconstructing from the example:
  E = 3.200 (a=0.100), Var = 7.716 (c=0.158) reproduce (v, d).
- Conservative (Theorem 5.2.2, p. 95): example with f_min = 21, a_min = -1.926 < -3/2,
  so v = 4 and critical value = (3 + 6/21)(1 - 2/4) F_0.95(3,4) = 1.64*F ~ 10.81.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from chang import combine
from chang import stats_dist as sd

print("=== Welch full data (Kirkwood and Chang 1998, Appendix B p.99) ===")
# a=0.100, b=0, c=0.158 (consistent with the paper example d=3.049, v=42.34)
w = combine.welch_df(0.100, 0.0, 0.158)
print(f"E={w['E']:.4f} (paper 3.200)   Var={w['Var']:.4f} (paper 7.716)")
print(f"df (v) = {w['df']:.2f}   (paper 42.34)")
print(f"cfact (d) = {w['cfact']:.4f}   (paper 3.049)")
# algebraic identity: v in terms of K=Var/E^2 must give (12K+2)/(3K-2)
K = w['Var'] / w['E'] ** 2
v_K = (12 * K + 2) / (3 * K - 2)
print(f"check of v via K=Var/E^2: {v_K:.2f} (must equal df)")

print("\n=== Welch critical value (full data) ===")
wc = combine.welch_critical_value(0.100, 0.0, 0.158, alpha=0.95)
print(f"crit = c*F_0.95(3, {wc['df']:.1f}) = {wc['crit']:.3f}")

print("\n=== Conservative (Theorem 5.2.2): reproduction of the example v=4 ===")
# The paper reports a_min=-1.926, f_min=21 => v=4, crit ~ 10.81. The final formula is
# checked directly with those intermediate values of the paper:
f_min = 21
v = 4.0
factor = (3.0 + 6.0 / f_min) * (1.0 - 2.0 / v)
crit_paper_path = factor * sd.f_crit(0.95, 3, v)
print(f"factor (3+6/21)(1-2/4) = {factor:.4f}  (paper 1.64)")
print(f"F_0.95(3,4) = {sd.f_crit(0.95, 3, 4):.4f}")
print(f"crit = {crit_paper_path:.3f}  (paper ~10.81)")

print("\n=== conservative_critical_value with data that force a_min < -3/2 ===")
# nf_pairs built with f_min=21 and large r so that z<0 strongly
# (two fits with many sections: f_i=21, large n_i, large r_i)
res = combine.conservative_critical_value([(120, 21), (120, 21)], alpha=0.95)
print(f"n={res['n']}, r={res['r']}, f_min={res['f_min']}, z={res['z']:.2f}, "
      f"a_min={res['a_min']:.3f}")
print(f"v={res['v']:.2f}, factor={res['factor']:.4f}, crit={res['crit']:.3f}")
print("(if a_min<-1.5 -> v=4; crit uses the conservative bound of Theorem 5.2.2)")

print("\n=== conservative with moderate kappas/df (typical case, a_min>-3/2) ===")
res2 = combine.conservative_critical_value([(40, 25), (60, 45)], alpha=0.95)
print(f"n={res2['n']}, r={res2['r']}, f_min={res2['f_min']}, a_min={res2['a_min']:.3f}, "
      f"K={res2['K']:.4f}, v={res2['v']:.2f}, crit={res2['crit']:.3f}")
print("comparison: 3*F_0.95(3, f_min) (equal kappas) =",
      round(3 * sd.f_crit(0.95, 3, res2['f_min']), 3),
      "-> the conservative value must be somewhat larger")
