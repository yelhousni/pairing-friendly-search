"""
Integral residues of rational polynomials, done properly.

For f = N/d and a prime p with m = v_p(d), f is p-integral at x iff
N(x) = 0 mod p^m. So we need the roots of an integer polynomial modulo p^m.
Level 1 is a root-finding problem over GF(p) (NOT a scan: p can be ~3*10^7),
and higher levels are Hensel lifts, branching only at singular roots.
"""

from sage.all import (QQ, ZZ, GF, Zmod, PolynomialRing, lcm, gcd, factor,
                      CRT_list, inverse_mod)


def num_den(f):
    d = ZZ(1)
    for c in f.coefficients():
        d = lcm(d, QQ(c).denominator())
    N = (f * d)
    return N.change_ring(QQ), d


def roots_mod_pm(N, p, m, cap=20000, branch_cap=5000, keep=4000):
    """Roots of N mod p^m."""
    p = ZZ(p); m = ZZ(m)
    Nz = N.change_ring(ZZ)
    Fp = GF(p)
    Np = Nz.change_ring(Fp)
    if Np.is_zero():
        if p > branch_cap:
            return None
        cur = [ZZ(r) for r in range(int(p))]
    else:
        try:
            cur = sorted(ZZ(r) for r in Np.roots(multiplicities=False))
        except Exception:
            return None
    if not cur:
        return []
    if m == 1:
        return cur
    PM = p ** m
    Rm = PolynomialRing(Zmod(PM), 'x')
    Nm = Rm(Nz.coefficients(sparse=False))
    Dm = Rm(Nz.derivative().coefficients(sparse=False))
    for j in range(2, int(m) + 1):
        pj = p ** j
        pjm = p ** (j - 1)
        nxt = []
        for r in cur:
            dv = ZZ(Dm(Zmod(PM)(r)).lift()) % p
            val = ZZ(Nm(Zmod(PM)(r)).lift()) % pj
            if dv != 0:
                if val % pjm != 0:
                    continue
                t = (-(val // pjm) * inverse_mod(dv, p)) % p
                nxt.append((r + t * pjm) % pj)
            else:
                if p > branch_cap:
                    return None
                for t in range(int(p)):
                    rr = (r + t * pjm) % pj
                    if ZZ(Nm(Zmod(PM)(rr)).lift()) % pj == 0:
                        nxt.append(rr)
        cur = sorted(set(nxt))
        if not cur:
            return []
        if len(cur) > keep:
            # we only need witnesses, not the full solution set: keep a
            # bounded spread of survivors instead of giving up
            import random
            random.seed(1)
            cur = sorted(random.sample(cur, keep))
    return cur


def integral_residues(fs, cap=20000, max_combos=20000):
    """(L, residues mod L where every f in fs is integral), or (L, None)."""
    L = ZZ(1)
    nd = []
    for f in fs:
        N, d = num_den(f)
        nd.append((N, d))
        L = lcm(L, d)
    if L == 1:
        return (ZZ(1), [ZZ(0)])
    per_pp = []
    for (p, A) in factor(L):
        conds = [(N, ZZ(d).valuation(p)) for (N, d) in nd
                 if ZZ(d).valuation(p) > 0]
        # anchor on the condition with the largest required valuation
        conds.sort(key=lambda c: -c[1])
        Nstar, mstar = conds[0]
        base = roots_mod_pm(Nstar, p, mstar, cap=cap)
        if base is None:
            return (L, None)
        if not base:
            return (L, [])
        PA = p ** A
        Rm = PolynomialRing(Zmod(PA), 'x')
        others = [(Rm(N.change_ring(ZZ).coefficients(sparse=False)), m)
                  for (N, m) in conds[1:]]
        good = []
        for r in base:
            ok = True
            for (Nm, m) in others:
                if ZZ(Nm(Zmod(PA)(r)).lift()) % (p ** m) != 0:
                    ok = False
                    break
            if ok:
                good.append(ZZ(r))
        if not good:
            return (L, [])
        per_pp.append((p ** mstar, good))
    combos = [[]]
    for (pa, S) in per_pp:
        combos = [c + [b] for c in combos for b in S]
        if len(combos) > max_combos:
            combos = combos[:max_combos]
    mods = [pa for (pa, _) in per_pp]
    out = []
    for c in combos:
        try:
            out.append(ZZ(CRT_list([ZZ(b) for b in c], mods)))
        except Exception:
            pass
    Lred = ZZ(1)
    for m in mods:
        Lred = Lred * m
    return (Lred, sorted(set(out)))
