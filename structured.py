"""
Structured Brezing-Weng search: w = zeta_ell^e * pi.

Empirical observation (verified against KSS16/18/32/36 and Scott-Guillevic 54):
every known KSS-type family arises from a primitive element of the very
special shape

    w = zeta_ell * pi,     pi in O_F,  F = Q(sqrt(-D)) a quadratic subfield
                                           of Q(zeta_ell)

For these, w^{n/2} lies in F (whenever zeta_ell^{n/2} is a root of unity of F),
so the minimal polynomial is a quadratic in x^{n/2}:

    r(x) = x^n - Tr(w^{n/2}) x^{n/2} + N(pi)^{n/2}

which explains the 625=5^4, 13^8, 343=7^3, 7^6, 3^9 constant terms of the
published families. This module sweeps (ell, D, pi, e) directly.
"""

from sage.all import (CyclotomicField, euler_phi, QQ, ZZ, matrix, vector,
                      PolynomialRing, gcd, lcm, cyclotomic_polynomial,
                      is_squarefree, divisors)

Rx = PolynomialRing(QQ, 'x')
X = Rx.gen()


def available_D(K, Dmax):
    out = {}
    for D in range(1, Dmax + 1):
        if not is_squarefree(D):
            continue
        e = K(-D)
        try:
            if e.is_square():
                out[D] = e.sqrt()
        except Exception:
            pass
    return out


def omega_of(D, s):
    """Generator of the ring of integers of Q(sqrt(-D)), given s = sqrt(-D)."""
    if D % 4 == 3:          # -D = 1 mod 4
        return (1 + s) / 2
    return s


def power_basis_matrix(w, n):
    cols = []; p = w.parent().one()
    for _ in range(n):
        cols.append(p.list()); p = p * w
    return matrix(QQ, cols).transpose()


def represents_primes(f, Lcap=200000, nvals=6):
    if f.degree() < 1 or f.leading_coefficient() <= 0:
        return False
    L = ZZ(1)
    for c in f.coefficients():
        L = lcm(L, QQ(c).denominator())
    if L > Lcap:
        return False
    res = [b for b in range(int(L)) if f(QQ(b)).denominator() == 1]
    if not res:
        return False
    g = ZZ(0)
    for b in res:
        for m in range(nvals):
            g = gcd(g, ZZ(f(QQ(b + m * L))))
            if g == 1:
                return True
    return g == 1


def valid(r, t, y, q, k, D):
    """Cheap in-loop filter. The primality/integrality audit is deferred to
    rigorous_check, because q can have a very large denominator (KSS32's q has
    denominator 2970292) and scanning its residues is far too slow here."""
    if q.leading_coefficient() <= 0:              return False
    if not q.is_irreducible():                    return False
    if (cyclotomic_polynomial(k)(t - 1) % r) != 0: return False
    if ((q + 1 - t) % r) != 0:                    return False
    if (4 * q - t ** 2) != D * y ** 2:            return False
    return True


def integral_residues_crt(f, pp_cap=2 * 10 ** 5):
    """Residues mod L (L = lcm of coefficient denominators) where f is
    integral, computed prime-power-wise and recombined by CRT. Avoids
    scanning L values when L is large but smooth."""
    from sage.all import factor, CRT_list
    L = ZZ(1)
    for c in f.coefficients():
        L = lcm(L, QQ(c).denominator())
    if L == 1:
        return (ZZ(1), [ZZ(0)])
    per_pp = []
    for (p, a) in factor(L):
        pa = p ** a
        if pa > pp_cap:
            return (L, None)              # give up, caller falls back
        S = [b for b in range(int(pa))
             if QQ(f(QQ(b))).denominator() % p != 0]
        if not S:
            return (L, [])                # never integral at integers
        per_pp.append((ZZ(pa), S))
    # recombine (cap the number of combinations)
    combos = [[]]
    for (pa, S) in per_pp:
        combos = [c + [b] for c in combos for b in S[:40]][:4000]
    mods = [pa for (pa, _) in per_pp]
    out = []
    for c in combos[:4000]:
        try:
            out.append(ZZ(CRT_list([ZZ(b) for b in c], mods)))
        except Exception:
            pass
    return (L, out)


def rigorous_check(r, t, y, q, k, D, nvals=8):
    """Full Freeman-Scott-Teske audit of a candidate family."""
    errs = []
    if not r.is_irreducible():                     errs.append('r reducible')
    if not q.is_irreducible():                     errs.append('q reducible')
    if q.leading_coefficient() <= 0:               errs.append('q lead<=0')
    if (cyclotomic_polynomial(k)(t - 1) % r) != 0: errs.append('r!|Phi_k(t-1)')
    if ((q + 1 - t) % r) != 0:                     errs.append('r!|q+1-t')
    if (4 * q - t ** 2) != D * y ** 2:             errs.append('CM fails')
    # integrality + Bunyakovsky content on a common progression
    M = ZZ(4)
    for f in (r, t, y, q):
        for c in f.coefficients():
            M = lcm(M, QQ(c).denominator())
    L, res = integral_residues_crt(q)
    if res is None:
        errs.append('integrality undecided (L=%s)' % L)
        return errs
    good = []
    for b in res:
        if all(QQ(f(QQ(b))).denominator() == 1 for f in (r, t, q)):
            good.append(b)
    if not good:
        errs.append('no integer x0 with r,t,q all integral')
        return errs
    gq = ZZ(0); gr = ZZ(0)
    for b in good[:20]:
        for m in range(nvals):
            x0 = QQ(b + m * L)
            if QQ(q(x0)).denominator() == 1:
                gq = gcd(gq, ZZ(q(x0)))
            if QQ(r(x0)).denominator() == 1:
                gr = gcd(gr, ZZ(r(x0)))
    if gq != 1: errs.append('q content %s' % gq)
    if gr != 1: errs.append('r content %s' % gr)
    return errs


def prim_int(f):
    L = ZZ(1)
    for c in f.coefficients():
        L = lcm(L, QQ(c).denominator())
    g = L * f
    cont = ZZ(0)
    for c in g.coefficients():
        cont = gcd(cont, ZZ(c))
    if cont > 1:
        g = g / cont
    if g.leading_coefficient() < 0:
        g = -g
    return g


def search_ell(ell, Dmax=40, B=8, exps=(1,), k_list=None, rho_cap=QQ(2),
               verbose=False):
    """All families from w = zeta_ell^e * (a + b*omega_D)."""
    K = CyclotomicField(ell); n = euler_phi(ell); z = K.gen()
    if n < 2:
        return []
    Ds = available_D(K, Dmax)
    if not Ds:
        return []
    ks = k_list if k_list else [d for d in divisors(ell) if d >= 3]
    out = []
    for D, s in Ds.items():
        om = omega_of(D, s)
        for a in range(-B, B + 1):
            for b in range(-B, B + 1):
                # b=0,a=1 is pi=1, i.e. w=zeta_ell: the CLASSICAL cyclotomic
                # Brezing-Weng/BLS baseline, which must be included.
                # conj(zeta*pi)=zeta^-1*conj(pi), so with e=1 both signs of b
                # are needed.
                if gcd(a, b) != 1:
                    continue
                if b == 0 and a != 1:
                    continue
                pi = a + b * om
                for e in exps:
                    w = z ** e * pi
                    M = power_basis_matrix(w, n)
                    try:
                        Minv = M.inverse()
                    except Exception:
                        continue
                    r = w.minpoly().change_ring(QQ)
                    if r.degree() != n:
                        continue
                    rI = prim_int(r)
                    for k in ks:
                        if ell % k != 0:
                            continue
                        step = ell // k
                        for mm in range(1, k + 1):
                            if gcd(mm, k) != 1:
                                continue
                            zk = z ** (step * mm)
                            t = Rx(list(Minv * vector(QQ, zk.list()))) + 1
                            yv = (zk - 1) / s
                            y = Rx(list(Minv * vector(QQ, yv.list())))
                            q = (t ** 2 + D * y ** 2) / 4
                            if q.degree() < 1:
                                continue
                            rho = QQ(q.degree()) / QQ(n)
                            if rho > rho_cap:
                                continue
                            if not valid(r, t, y, q, k, D):
                                continue
                            out.append(dict(k=int(k), ell=int(ell), D=int(D),
                                            n=int(n), a=int(a), b=int(b),
                                            e=int(e), rho=rho, r=r, rI=rI,
                                            t=t, y=y, q=q,
                                            Npi=ZZ(pi.norm()) if hasattr(pi, 'norm') else None))
    return out
