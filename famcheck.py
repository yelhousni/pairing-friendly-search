"""
Test whether a potential family is an actual family, in the sense of
Gasnier-Guillevic Definition 2.2/2.4: Q irreducible with positive leading
coefficient, integer-valued somewhere, and with gcd of its values equal to 1.

This is the Bunyakovsky criterion, NOT an explicit prime seed. For deg R ~ 56
and deg Q ~ 60 the density of simultaneous primes is ~1/35000, so demanding a
seed is infeasible and is not what the definition requires. Use verify.py when
an explicit curve is actually wanted.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sage.all import QQ, ZZ, gcd, factor
from subfield import build, audit
from pplift import integral_residues


def family_check(k, D, a, b, nvals=4, max_residues=60):
    """Returns (status, info). status == 'FAMILY' iff it is an actual family."""
    f = build(k, D, a, b)
    if f is None:
        return ('build-none', None)
    errs = audit(f)
    if errs:
        return ('audit:' + ','.join(errs), None)
    R, T, Q = f['R'], f['T'], f['Q']
    L, res = integral_residues([R, T, Q])
    if res is None:
        return ('UNDECIDED', None)
    if not res:
        return ('never-integral', None)
    for b0 in res[:max_residues]:
        g = ZZ(0)
        for i in range(nvals):
            g = gcd(g, ZZ(Q(QQ(b0 + i * L))))
            if g == 1:
                break
        if g == 1:
            gr = ZZ(0)
            for i in range(nvals):
                gr = gcd(gr, ZZ(R(QQ(b0 + i * L))))
            return ('FAMILY', dict(rho=f['rho'], L=L, b0=b0,
                                   degR=R.degree(), degQ=Q.degree(),
                                   contentR=gr, R=R, T=T, Q=Q))
    return ('content(Q)>1 on %d residues' % len(res), None)


if __name__ == '__main__':
    k, D, a, b = (int(v) for v in sys.argv[1:5])
    st, info = family_check(k, D, a, b)
    print(f"k={k} D={D} alpha={a}{b:+d}w -> {st}")
    if info:
        print(f"  rho={info['rho']}  deg R={info['degR']}  deg Q={info['degQ']}")
        print(f"  Q integral and content-free on x = {info['L']}*u + {info['b0']}")
        print(f"  content(R) = {info['contentR']} (divide R by this)")
