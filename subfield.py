"""
The subfield method of Gasnier-Guillevic (SIAM J. Appl. Alg. Geom. 9(2), 2025),
implemented in closed form so it can be pushed past their k <= 54 range.

Setup: k >= 1, D squarefree > 0, F = Q(sqrt(-D)), C_k = Q(zeta_k),
K = F.C_k. Take theta = alpha * zeta_k with alpha in O_F, theta primitive in K.
Let e = minimal divisor of k with zeta_k^e in F.

When sqrt(-D) is NOT in C_k, F cap C_k = Q, so the only roots of unity in F
are +-1 and e = k/2 (k even) or e = k (k odd); then [K:Q] = 2*phi(k).

Let beta = theta^e in F. If beta generates F then {1, beta} is a Q-basis of F,
so there are DEGREE <= 1 rational polynomials P1,P2,P3 with
    P1(beta) = 1/alpha,  P2(beta) = 1/(alpha*sqrt(-D)),  P3(beta) = 1/sqrt(-D)
and then (Gasnier-Guillevic Thm 3.1)
    T(X) = P1(X^e)*X + 1                 -> theta maps to zeta_k + 1
    Y(X) = P2(X^e)*X - P3(X^e)           -> theta maps to (zeta_k-1)/sqrt(-D)
    Q(X) = (T^2 + D*Y^2)/4,   rho = deg Q / deg R = (2e+2)/(2 phi(k))
provided e+1 < deg R, in which case no reduction mod R is needed.

R is obtained in closed form: the conjugates of theta are alpha*zeta_k^j and
conj(alpha)*zeta_k^j over j coprime to k, so
    R(X) = [alpha^phi * Phi_k(X/alpha)] * [conj(alpha)^phi * Phi_k(X/conj(alpha))]
which is a monic integer polynomial of degree 2*phi(k).
"""

from sage.all import (QQ, ZZ, PolynomialRing, QuadraticField, cyclotomic_polynomial,
                      euler_phi, gcd, lcm, is_squarefree, CyclotomicField, matrix,
                      vector, divisors)

Rx = PolynomialRing(QQ, 'x')
X = Rx.gen()


def sqrt_minus_D_in_Ck(k, D):
    """Is sqrt(-D) in Q(zeta_k)?"""
    L = CyclotomicField(k)
    try:
        return L(-D).is_square()
    except Exception:
        return False


def omega_gen(F, D):
    """Generator of O_F: (1+sqrt(-D))/2 if -D = 1 mod 4, else sqrt(-D)."""
    s = F.gen()                       # sqrt(-D)
    if D % 4 == 3:
        return (1 + s) / 2
    return s


def build(k, D, a, b, check_primitive=True):
    """Build the potential family for theta = (a + b*omega) * zeta_k."""
    if D % 4 == 3:
        pass
    F = QuadraticField(-D, 'sq')
    s = F.gen()
    om = omega_gen(F, D)
    alpha = a + b * om
    if alpha == 0:
        return None
    phi = euler_phi(k)

    # e = minimal divisor of k with zeta_k^e in F. With sqrt(-D) not in C_k the
    # only roots of unity in F are +-1.
    if k % 2 == 0:
        e = k // 2
        zke = -1
    else:
        e = k
        zke = 1

    beta = (alpha ** e) * zke
    if beta in QQ:
        return None                    # theta^e not primitive in F

    # coordinates in the Q-basis {1, beta} of F
    def coords(v):
        v = F(v)
        M = matrix(QQ, [list(F(1)), list(beta)]).transpose()
        try:
            c = M.solve_right(vector(QQ, list(v)))
        except Exception:
            return None
        return (QQ(c[0]), QQ(c[1]))

    c1 = coords(1 / alpha)
    c2 = coords(1 / (alpha * s))
    c3 = coords(1 / s)
    if c1 is None or c2 is None or c3 is None:
        return None
    P1 = Rx([c1[0], c1[1]])
    P2 = Rx([c2[0], c2[1]])
    P3 = Rx([c3[0], c3[1]])

    T = P1(X ** e) * X + 1
    Y = P2(X ** e) * X - P3(X ** e)
    Q = (T ** 2 + D * Y ** 2) / 4

    # R = [alpha^phi Phi_k(X/alpha)] * [conj(alpha)^phi Phi_k(X/conj(alpha))]
    Phi = cyclotomic_polynomial(k)
    co = Phi.coefficients(sparse=False)
    FX = PolynomialRing(F, 'x'); XF = FX.gen()
    ab = alpha.conjugate()
    A = sum(F(co[i]) * alpha ** (phi - i) * XF ** i for i in range(phi + 1))
    B = sum(F(co[i]) * ab ** (phi - i) * XF ** i for i in range(phi + 1))
    RR = A * B
    coeffs = []
    for c in RR.coefficients(sparse=False):
        cc = F(c)
        l = list(cc)
        if l[1] != 0:
            return None                # should be rational
        coeffs.append(QQ(l[0]))
    R = Rx(coeffs)

    if check_primitive and not R.is_irreducible():
        return None
    if R.degree() != 2 * phi:
        return None
    if e + 1 >= R.degree():
        return None                    # reduction needed; outside Thm 3.1 regime
    rho = QQ(Q.degree()) / QQ(R.degree())
    return dict(k=k, D=D, a=a, b=b, e=e, R=R, T=T, Y=Y, Q=Q, rho=rho,
                alpha=alpha, phi=phi)


def audit(f):
    """Symbolic conditions of Gasnier-Guillevic Definition 2.4."""
    R, T, Y, Q, k, D = f['R'], f['T'], f['Y'], f['Q'], f['k'], f['D']
    errs = []
    if not R.is_irreducible():                      errs.append('R reducible')
    if (cyclotomic_polynomial(k)(T - 1) % R) != 0:  errs.append('R !| Phi_k(T-1)')
    if ((Q + 1 - T) % R) != 0:                      errs.append('R !| Q+1-T')
    if (4 * Q - T ** 2) != D * Y ** 2:              errs.append('CM fails')
    if not Q.is_irreducible():                      errs.append('Q reducible')
    if Q.leading_coefficient() <= 0:                errs.append('Q lead<=0')
    return errs


def rho_formula(k, D):
    """Thm 3.3 cases 4,5: sqrt(-D) not in C_k."""
    phi = euler_phi(k)
    if k % 2 == 0:
        return QQ(k // 2 + 1) / QQ(phi)
    return QQ(k + 1) / QQ(phi)
