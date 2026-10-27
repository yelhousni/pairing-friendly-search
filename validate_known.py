"""
Validation suite: reproduce published families from scratch.

Nothing in this repository should be trusted until this passes. Every claim
the search makes rests on the engine agreeing with the literature on families
whose answers are already known.

Run:  sage validate_known.py          (a few minutes)
      sage validate_known.py quick    (skip the slow k=48 case)
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sage.all import (QQ, ZZ, PolynomialRing, CyclotomicField, euler_phi,
                      cyclotomic_polynomial, lcm, gcd)
from subfield import build as sub_build, audit as sub_audit
from structured import search_ell
from normalize import verify_family, progressions, _denom
from fast import sweep_fast
from l3 import l3_families

Rx = PolynomialRing(QQ, 'x'); x = Rx.gen()
QUICK = len(sys.argv) > 1 and sys.argv[1] == 'quick'
results = []


def report(name, ok, detail=''):
    results.append((name, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ''),
          flush=True)


# ---------------------------------------------------------------- 1. BN
print("\n1. Barreto-Naehrig (k=12, D=3, rho=1)", flush=True)
t0 = time.time()
fams = l3_families(12, 6 * x ** 2, [12], Dmax=6)
hit = [f for f in fams if f['rho'] == 1]
ok = any(str(f['r']) == '36*x^4 + 36*x^3 + 18*x^2 + 6*x + 1' and
         str(f['t']) == '6*x^2 + 1' and
         str(f['q']) == '36*x^4 + 36*x^3 + 24*x^2 + 6*x + 1' for f in hit)
report('BN from Phi_12(6x^2) factorisation', ok, f'[{time.time()-t0:.1f}s]')

# BN is also reachable as a primitive element of Q(zeta_12)
t0 = time.time()
res = sweep_fast(12, 12, Dmax=6, specs=[('full', 3, None)])
ok = any(v[0][0] == 1 for v in res.values())
report('BN-equivalent rho=1 found by dense sweep of Q(zeta_12)', ok,
       f'[{time.time()-t0:.1f}s]')


# ------------------------------------------- 2. cyclotomic / KSS families
print("\n2. Cyclotomic and KSS families (structured search w = zeta_ell * pi)",
      flush=True)
targets = [(16, 16, QQ(5) / 4, 'KSS16'), (18, 18, QQ(4) / 3, 'KSS18'),
           (36, 36, QQ(7) / 6, 'KSS36'), (54, 54, QQ(10) / 9, 'SG54'),
           (24, 8, QQ(5) / 4, 'KSS8'), (24, 24, QQ(5) / 4, 'BLS24'),
           (12, 12, QQ(3) / 2, 'BLS12')]
if not QUICK:
    targets.append((48, 48, QQ(9) / 8, 'BLS48'))
for ell, k, want, label in targets:
    t0 = time.time()
    cands = [f for f in search_ell(ell, Dmax=20, B=6, k_list=[k])
             if _denom(f['q']) <= 200000]
    cands.sort(key=lambda f: (f['rho'], abs(f['a']) + abs(f['b'])))
    found = None
    for f in cands[:400]:
        if not progressions(f['rI'], f['t'], f['q']):
            continue
        if verify_family(f['rI'], f['t'], f['q'], f['D'], k, umax=1500, want=1):
            found = f
            break
    ok = found is not None and found['rho'] == want
    report(f'{label}: rho = {want}', ok,
           f"got {found['rho'] if found else None}  [{time.time()-t0:.1f}s]")


# --------------------------------------------- 3. FST Construction 6.7
print("\n3. Freeman-Scott-Teske Construction 6.7 at k=72 (D=2)", flush=True)
t0 = time.time()
k = 72; ell = lcm(8, k); a = ell // k
t_f = x ** a + 1
y_f = (1 - x ** a) * (x ** (5 * ell // 24) + x ** (ell // 8) - x ** (ell // 24)) / 2
q_f = (2 * (x ** a + 1) ** 2 +
       (1 - x ** a) ** 2 * (x ** (5 * ell // 24) + x ** (ell // 8) - x ** (ell // 24)) ** 2) / 8
K = CyclotomicField(72); z = K.gen(); s = K(-2).sqrt()
t_m = Rx(list((z + 1).list()))
y_m = Rx(list(((z - 1) / s).list()))
q_m = (t_m ** 2 + 2 * y_m ** 2) / 4
report('closed form 6.7 == direct Brezing-Weng with D=2',
       t_f == t_m and y_f == y_m and q_f == q_m)
rho = QQ(q_f.degree()) / QQ(euler_phi(k))
report('rho = 4/3 matches the (5k/12+2)/phi(k) formula',
       rho == QQ(4) / 3 and rho == QQ(5 * k // 12 + 2) / QQ(euler_phi(k)),
       f'rho={rho}  [{time.time()-t0:.1f}s]')


# ------------------------------- 4. Gasnier-Guillevic subfield method
print("\n4. Gasnier-Guillevic subfield method, Example 3.4 (k=22, D=7)",
      flush=True)
t0 = time.time()
f = sub_build(22, 7, 0, 1)          # alpha = omega, norm 2
T_p = (x ** 12 + 45 * x + 46) / 46
Y_p = (x ** 12 - 4 * x ** 11 - 47 * x - 134) / 322
R_p = (x ** 20 - x ** 19 - x ** 18 + 3 * x ** 17 - x ** 16 - 5 * x ** 15 +
       7 * x ** 14 + 3 * x ** 13 - 17 * x ** 12 + 11 * x ** 11 + 23 * x ** 10 +
       22 * x ** 9 - 68 * x ** 8 + 24 * x ** 7 + 112 * x ** 6 - 160 * x ** 5 -
       64 * x ** 4 + 384 * x ** 3 - 256 * x ** 2 - 512 * x + 1024) / 23
Q_p = (x ** 24 - x ** 23 + 2 * x ** 22 + 67 * x ** 13 + 94 * x ** 12 +
       134 * x ** 11 + 2048 * x ** 2 + 5197 * x + 4096) / 7406
ok = (f is not None and not sub_audit(f) and f['rho'] == QQ(6) / 5 and
      f['T'] == T_p and f['Q'] == Q_p and f['R'] == 23 * R_p and
      (f['Y'] == Y_p or f['Y'] == -Y_p))
report('R, T, Y, Q all match the published polynomials', ok,
       f"rho={f['rho'] if f else None}  [{time.time()-t0:.1f}s]")


# ----------------------------------------------------------------- done
print("\n" + "=" * 62)
npass = sum(1 for _, ok in results if ok)
print(f"{npass}/{len(results)} checks passed")
print("=" * 62)
sys.exit(0 if npass == len(results) else 1)
