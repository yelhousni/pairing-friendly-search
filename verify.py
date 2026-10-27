"""
Gold-standard verification of a candidate family (r,t,q,D,k).

Beyond the symbolic conditions we look for actual integer seeds x0 giving
primes, and then check the arithmetic facts directly:
  - q = q(x0) is prime, r = r(x0) is prime
  - r | q+1-t   (so a curve of trace t over F_q has an order-r subgroup)
  - t^2-4q = -D*v^2  (CM equation solvable => the curve exists, via CM method)
  - |t| <= 2*sqrt(q) (Hasse)
  - the multiplicative order of q modulo r is exactly k (embedding degree)
"""
from sage.all import (QQ, ZZ, is_prime, Mod, gcd, sqrt, divisors)


def check_seed(r, t, q, D, k, x0):
    x0 = QQ(x0)
    rv, tv, qv = r(x0), t(x0), q(x0)
    if any(QQ(v).denominator() != 1 for v in (rv, tv, qv)):
        return None
    rv, tv, qv = ZZ(rv), ZZ(tv), ZZ(qv)
    if rv <= 1 or qv <= 3:
        return None
    if not is_prime(qv) or not is_prime(rv):
        return None
    if (qv + 1 - tv) % rv != 0:
        return None
    disc = tv ** 2 - 4 * qv
    if disc >= 0:
        return None
    if (-disc) % D != 0:
        return None
    v2 = (-disc) // D
    if not ZZ(v2).is_square():
        return None
    if tv ** 2 > 4 * qv:
        return None
    # embedding degree: order of q mod r must be exactly k
    if gcd(qv, rv) != 1:
        return None
    if Mod(qv, rv).multiplicative_order() != k:
        return None
    return dict(x0=x0, r=rv, t=tv, q=qv, h=(qv + 1 - tv) // rv,
                rbits=rv.nbits(), qbits=qv.nbits())


def find_seeds(r, t, q, D, k, xmax=20000, want=3, denom_hint=None):
    """Scan integer (and small-denominator rational) seeds."""
    from sage.all import lcm
    L = ZZ(1)
    for f in (r, t, q):
        for c in f.coefficients():
            L = lcm(L, QQ(c).denominator())
    out = []
    for x0 in range(1, xmax):
        got = check_seed(r, t, q, D, k, x0)
        if got:
            out.append(got)
            if len(out) >= want:
                return out
        got = check_seed(r, t, q, D, k, -x0)
        if got:
            out.append(got)
            if len(out) >= want:
                return out
    return out


def rho_empirical(seed):
    return float(seed['q'].nbits()) / float(seed['r'].nbits())
