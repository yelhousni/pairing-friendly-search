"""
Family normalization + gold-standard verification.

A Brezing-Weng family is only defined up to the reparametrization
x -> L*u + b, and up to dividing r by a constant. Both are essential:
the standard KSS18 family is

    r(x) = (x^6+37x^3+343)/343,   x = 14 mod 42

i.e. the raw minimal polynomial x^6+37x^3+343 has content 343 on the relevant
progression and must be divided by it. (q may NOT be divided by a content: a
common factor in q's values means q is never prime, so such a progression is
simply rejected.)

Final arbiter is not symbolic at all: we look for integer seeds u giving
r, q both prime with r | q+1-t, the CM equation solvable, and the
multiplicative order of q mod r exactly k.
"""

from sage.all import (QQ, ZZ, lcm, gcd, factor, CRT_list, PolynomialRing,
                      is_prime, Mod)

Ru = PolynomialRing(QQ, 'u')
U = Ru.gen()


def _denom(f):
    L = ZZ(1)
    for c in f.coefficients():
        L = lcm(L, QQ(c).denominator())
    return L


def _content(f):
    g = ZZ(0)
    for c in f.coefficients():
        g = gcd(g, ZZ(c))
    return g if g != 0 else ZZ(1)


def progressions_crt(r, t, q, pp_cap=20000, max_combos=400):
    """Same as progressions() but finds the admissible residues prime-power
    wise and recombines by CRT, instead of scanning all L values. Needed once
    the denominators get large (subfield-method Q can have huge denominators)."""
    from sage.all import factor, CRT_list
    L = ZZ(1)
    for f in (r, t, q):
        L = lcm(L, _denom(f))
    if L == 1:
        residues = [ZZ(0)]
    else:
        per_pp = []
        for (p, a) in factor(L):
            pa = p ** a
            if pa > pp_cap:
                return []
            S = [b for b in range(int(pa))
                 if all(QQ(f(QQ(b))).denominator() % p != 0 for f in (r, t, q))]
            if not S:
                return []
            per_pp.append((ZZ(pa), S))
        combos = [[]]
        for (pa, S) in per_pp:
            combos = [c + [b] for c in combos for b in S][:max_combos]
        mods = [pa for (pa, _) in per_pp]
        residues = []
        for c in combos:
            try:
                residues.append(ZZ(CRT_list([ZZ(b) for b in c], mods)))
            except Exception:
                pass
    out = []
    seen = set()
    for b in residues:
        sub = L * U + b
        rb, tb, qb = r(sub), t(sub), q(sub)
        if any(_denom(f) != 1 for f in (rb, tb, qb)):
            continue
        rb = Ru(rb); tb = Ru(tb); qb = Ru(qb)
        if _content(qb) != 1:
            continue
        cr = _content(rb)
        rn = rb / cr
        if rn.leading_coefficient() < 0:
            rn = -rn
        key = (str(rn), str(tb), str(qb))
        if key in seen:
            continue
        seen.add(key)
        out.append((L, ZZ(b), rn, tb, qb))
    return out


def progressions(r, t, q, extra=(1, 2, 3, 4, 6, 12), scan_cap=200000):
    """Yield (L, b, r_n, t_n, q_n) reparametrizations x = L*u+b in which
    r,t,q are integral, q is content-free, and r has been content-normalized."""
    base = ZZ(1)
    for f in (r, t, q):
        base = lcm(base, _denom(f))
    out = []
    seen = set()
    for mult in extra:
        L = base * mult
        if L > scan_cap:
            continue
        for b in range(int(L)):
            sub = L * U + b
            rb, tb, qb = r(sub), t(sub), q(sub)
            if any(_denom(f) != 1 for f in (rb, tb, qb)):
                continue
            rb = Ru(rb); tb = Ru(tb); qb = Ru(qb)
            if _content(qb) != 1:
                continue
            cr = _content(rb)
            rn = rb / cr
            if rn.leading_coefficient() < 0:
                rn, = (-rn,)
            key = (str(rn), str(tb), str(qb))
            if key in seen:
                continue
            seen.add(key)
            out.append((L, ZZ(b), rn, tb, qb))
        if out:
            break
    return out


def check_seed(rn, tn, qn, D, k, u0):
    rv, tv, qv = rn(QQ(u0)), tn(QQ(u0)), qn(QQ(u0))
    if any(QQ(v).denominator() != 1 for v in (rv, tv, qv)):
        return None
    rv, tv, qv = ZZ(rv), ZZ(tv), ZZ(qv)
    if rv <= 3 or qv <= 3:
        return None
    if not is_prime(qv) or not is_prime(rv):
        return None
    N = qv + 1 - tv                      # curve order
    if N % rv != 0:
        return None
    disc = tv ** 2 - 4 * qv
    if disc >= 0 or (-disc) % D != 0:
        return None
    if not ZZ((-disc) // D).is_square():
        return None
    if gcd(qv, rv) != 1:
        return None
    if Mod(qv, rv).multiplicative_order() != k:
        return None
    return dict(u0=ZZ(u0), r=rv, t=tv, q=qv, h=N // rv,
                rbits=rv.nbits(), qbits=qv.nbits(),
                rho=float(qv.nbits()) / float(rv.nbits()))


def verify_family(r, t, q, D, k, umax=4000, want=2, use_crt=False):
    """Try every admissible progression; return seeds that pass all checks."""
    best = []
    progs = progressions_crt(r, t, q) if use_crt else progressions(r, t, q)
    for (L, b, rn, tn, qn) in progs:
        seeds = []
        for u0 in range(-umax, umax):
            s = check_seed(rn, tn, qn, D, k, u0)
            if s:
                s['L'] = L; s['b'] = b
                s['r_poly'] = rn; s['t_poly'] = tn; s['q_poly'] = qn
                seeds.append(s)
                if len(seeds) >= want:
                    break
        if seeds:
            best = seeds
            break
    return best
