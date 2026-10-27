"""
Layer 3: non-cyclotomic Brezing-Weng.

Instead of K = Q(zeta_ell), take K = Q[x]/(r) where r is an irreducible factor
of Phi_ell(u(x)) for some polynomial u. Then u(x) mod r is automatically a
primitive ell-th root of unity, so K contains zeta_ell. This is the layer
containing the rho = 1 families: BN is ell=12, u(x)=6x^2, where
Phi_12(6x^2) = 1296x^8-36x^4+1 factors into two quartics, one of which is
r(x) = 36x^4+36x^3+18x^2+6x+1.

Normalization: x -> lambda*x+mu is a family-preserving change of variable, so
a quadratic u can be taken as a*x^2+c, a cubic as a*x^3+c*x+d, etc.
"""

from sage.all import (QQ, ZZ, PolynomialRing, NumberField, cyclotomic_polynomial,
                      euler_phi, gcd, lcm, is_squarefree, divisors)
import itertools

Rx = PolynomialRing(QQ, 'x')
X = Rx.gen()


def u_candidates(deg, arange, crange):
    """Monic-free u of given degree, depressed (no second-highest term)."""
    if deg == 2:
        for a in arange:
            if a == 0:
                continue
            for c in crange:
                yield Rx([c, 0, a])
    elif deg == 3:
        for a in arange:
            if a == 0:
                continue
            for c in crange:
                for d in crange:
                    yield Rx([d, c, 0, a])
    elif deg == 4:
        for a in arange:
            if a == 0:
                continue
            for c in crange:
                for d in crange:
                    for e in crange:
                        yield Rx([e, d, c, 0, a])


def prim_int(f):
    """Scale f by a rational constant to make it integral and primitive,
    with positive leading coefficient. Same roots, so the same number field."""
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


def l3_families(ell, u, k_list, Dmax=20, rho_cap=QQ(2), require_factor=True):
    """All families from irreducible factors of Phi_ell(u(x))."""
    P = cyclotomic_polynomial(ell)(u)
    if P.degree() < 2:
        return []
    fac = P.factor()
    if require_factor and len(fac) == 1 and fac[0][1] == 1:
        return []                      # no nontrivial factorization
    out = []
    phi = euler_phi(ell)
    for (rm, mult) in fac:
        d = rm.degree()
        if d < phi or d > 40:
            continue
        rm = rm.monic()
        try:
            K = NumberField(rm, 'a')
        except Exception:
            continue
        r = prim_int(rm)        # family's r: primitive integral, same roots
        a = K.gen()
        zl = u(a)                      # primitive ell-th root of unity in K
        # available D
        Ds = {}
        for D in range(1, Dmax + 1):
            if not is_squarefree(D):
                continue
            try:
                e = K(-D)
                if e.is_square():
                    Ds[D] = e.sqrt()
            except Exception:
                pass
        if not Ds:
            continue
        for k in k_list:
            if ell % k != 0:
                continue
            step = ell // k
            for m in range(1, k + 1):
                if gcd(m, k) != 1:
                    continue
                zk = zl ** (step * m)
                t = Rx(list((zk + 1).list()))
                for D, s in Ds.items():
                    y = Rx(list(((zk - 1) / s).list()))
                    q = (t ** 2 + D * y ** 2) / 4
                    if q.degree() < 1:
                        continue
                    rho = QQ(q.degree()) / QQ(d)
                    if rho > rho_cap:
                        continue
                    if q.leading_coefficient() <= 0 or not q.is_irreducible():
                        continue
                    if not represents_primes(q) or not represents_primes(r):
                        continue
                    # audit
                    if (cyclotomic_polynomial(k)(t - 1) % r) != 0:
                        continue
                    if ((q + 1 - t) % r) != 0:
                        continue
                    if (4 * q - t ** 2) != D * y ** 2:
                        continue
                    out.append(dict(k=int(k), ell=int(ell), D=int(D), d=int(d),
                                    u=str(u), rho=rho, r=r, t=t, y=y, q=q))
    return out
