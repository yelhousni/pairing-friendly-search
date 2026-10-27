"""
Fast generalized Brezing-Weng sweep (two-pass).

Of the validity conditions, three are automatic by construction:
  t = zeta_k+1 mod r  =>  Phi_k(t-1) = 0 mod r
  q = (t^2+Dy^2)/4 = zeta_k mod r  =>  q+1-t = 0 mod r
  r = minpoly(w)  =>  r irreducible
  4q-t^2 = D y^2  holds identically in Q[x]
So the only nontrivial filters are on q: irreducible, positive leading
coefficient, represents primes.

Pass A sweeps every (w, zeta_k, D) and records only rho = deg q / n (cheap).
Pass B validates candidates in increasing rho order until valid ones are found.
"""

from sage.all import (CyclotomicField, euler_phi, QQ, ZZ, matrix, vector,
                      PolynomialRing, gcd, lcm, cyclotomic_polynomial,
                      is_squarefree)
import itertools
from collections import defaultdict

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


def primitive_kth_roots(K, ell, k):
    z = K.gen(); step = ell // k
    return [z ** (step * m) for m in range(1, k + 1) if gcd(m, k) == 1]


def power_basis_matrix(w, n):
    cols = []; p = w.parent().one()
    for _ in range(n):
        cols.append(p.list()); p = p * w
    return matrix(QQ, cols).transpose()


def enumerate_w(n, coeff_range, max_terms, max_idx=None):
    """Sparse mode: at most max_terms nonzero coefficients."""
    top = n if max_idx is None else min(n, max_idx + 1)
    idxs = list(range(1, top))
    C = [c for c in range(-coeff_range, coeff_range + 1) if c != 0]
    for nterms in range(1, max_terms + 1):
        for positions in itertools.combinations(idxs, nterms):
            for vals in itertools.product(C, repeat=nterms):
                if vals[0] < 0:
                    continue
                g = 0
                for c in vals:
                    g = gcd(g, c)
                if g > 1:
                    continue
                yield positions, vals


def enumerate_w_full(n, coeff_range, max_idx=None):
    """Dense mode: every coefficient vector (c_1..c_{n-1}) in [-C,C].

    Every element of K is equivalent, under the family-preserving affine map
    w -> lambda*w+mu, to a primitive integer vector with zero constant term and
    first nonzero coefficient positive. So this enumerates ALL primitive
    elements of K up to equivalence, within the coefficient box.
    """
    top = n if max_idx is None else min(n, max_idx + 1)
    m = top - 1
    C = list(range(-coeff_range, coeff_range + 1))
    for vec in itertools.product(C, repeat=m):
        first = 0
        for c in vec:
            if c != 0:
                first = c
                break
        if first <= 0:
            continue                      # zero vector, or normalize sign
        g = 0
        for c in vec:
            g = gcd(g, c)
        if g > 1:
            continue                      # normalize scaling
        positions = tuple(i + 1 for i, c in enumerate(vec) if c != 0)
        vals = tuple(c for c in vec if c != 0)
        yield positions, vals


def integral_residues(f, Lcap=200000):
    """The integers where f is integral form a union of residue classes mod L,
    L = lcm of the denominators of the coefficients. Return (L, [residues])."""
    L = ZZ(1)
    for c in f.coefficients():
        L = lcm(L, QQ(c).denominator())
    if L > Lcap:
        return (L, [])
    res = [b for b in range(int(L)) if f(QQ(b)).denominator() == 1]
    return (L, res)


def represents_primes(f, nvals=6):
    """Bunyakovsky-style: irreducible (checked by caller), positive leading
    coefficient, integer-valued on some progression, and gcd of its values
    on that progression equal to 1."""
    if f.degree() < 1 or f.leading_coefficient() <= 0:
        return False
    L, res = integral_residues(f)
    if not res:
        return False
    g = ZZ(0)
    for b in res:
        for m in range(nvals):
            g = gcd(g, ZZ(f(QQ(b + m * L))))
            if g == 1:
                return True
    return g == 1


def fully_valid(r, t, y, q, k, D):
    errs = []
    if q.leading_coefficient() <= 0:      errs.append('q lead<=0')
    elif not q.is_irreducible():          errs.append('q reducible')
    elif not represents_primes(q):        errs.append('q !primes')
    if errs:
        return errs
    if not r.is_irreducible():            errs.append('r reducible')
    if (cyclotomic_polynomial(k)(t - 1) % r) != 0: errs.append('r!|Phi_k(t-1)')
    if ((q + 1 - t) % r) != 0:            errs.append('r!|q+1-t')
    if (4 * q - t ** 2) != D * y ** 2:    errs.append('CM fails')
    return errs


def sweep_fast(k, ell, coeff_range=2, max_terms=3, Dmax=20, max_idx=None,
               rho_cap=QQ(2), n_rho_levels=14, per_level=400, keep=4,
               specs=None, roots_limit=None):
    """specs: list of ('full', C, None) or ('sparse', C, max_terms) passes.
    If None, a single sparse pass with (coeff_range, max_terms)."""
    if ell % k != 0:
        return {}
    K = CyclotomicField(ell); n = euler_phi(ell)
    if n < 2:
        return {}
    Ds = available_D(K, Dmax)
    if not Ds:
        return {}
    roots = primitive_kth_roots(K, ell, k)
    if not roots:
        return {}
    if roots_limit is not None:
        roots = roots[:roots_limit]
    if specs is None:
        specs = [('sparse', coeff_range, max_terms)]

    def wgen():
        seen = set()
        for (mode, C, mt) in specs:
            it = (enumerate_w_full(n, C, max_idx) if mode == 'full'
                  else enumerate_w(n, C, mt, max_idx))
            for positions, vals in it:
                key = (positions, vals)
                if key in seen:
                    continue
                seen.add(key)
                yield positions, vals
    zvecs = [vector(QQ, z.list()) for z in roots]
    yvecs = {D: [vector(QQ, ((z - 1) / s).list()) for z in roots]
             for D, s in Ds.items()}
    z = K.gen()

    # ---- Pass A: cheap rho only
    buckets = {D: defaultdict(list) for D in Ds}   # rho -> [(positions,vals,ri)]
    for positions, vals in wgen():
        w = sum(ZZ(c) * z ** i for i, c in zip(positions, vals))
        M = power_basis_matrix(w, n)
        if M.determinant() == 0:
            continue
        Minv = M.inverse()
        for ri in range(len(roots)):
            t = Rx(list(Minv * zvecs[ri])) + 1
            for D in Ds:
                y = Rx(list(Minv * yvecs[D][ri]))
                q = (t ** 2 + D * y ** 2) / 4
                dq = q.degree()
                if dq < 1:
                    continue
                rho = QQ(dq) / QQ(n)
                if rho > rho_cap:
                    continue
                b = buckets[D]
                if rho in b:
                    if len(b[rho]) < per_level:
                        b[rho].append((positions, vals, ri))
                else:
                    b[rho] = [(positions, vals, ri)]
                    if len(b) > n_rho_levels:
                        del b[max(b)]

    # ---- Pass B: validate in increasing rho order
    out = {}
    for D, b in buckets.items():
        s = Ds[D]
        found = []
        for rho in sorted(b):
            for (positions, vals, ri) in b[rho]:
                w = sum(ZZ(c) * z ** i for i, c in zip(positions, vals))
                r = w.minpoly().change_ring(QQ)
                if r.degree() != n:
                    continue
                M = power_basis_matrix(w, n)
                Minv = M.inverse()
                zeta = roots[ri]
                t = Rx(list(Minv * vector(QQ, zeta.list()))) + 1
                y = Rx(list(Minv * vector(QQ, ((zeta - 1) / s).list())))
                q = (t ** 2 + D * y ** 2) / 4
                if fully_valid(r, t, y, q, k, D):
                    continue
                found.append((rho, r, t, y, q, (positions, vals, ri)))
                if len(found) >= keep:
                    break
            if found:
                break
        if found:
            out[D] = found
    return out
