# pairing-friendly-search

A SageMath toolkit for searching for families of pairing-friendly elliptic
curves, covering the Brezing–Weng / KSS / subfield constructions under one
roof, with a validation suite that reproduces the published families exactly.

Given an embedding degree `k`, a CM discriminant `D` and a search space of
number fields, it enumerates candidate families, filters them down to the ones
that are genuinely families (not merely *potential* families), and will
produce explicit curves of cryptographic size where they exist.

## Validation

Nothing should be trusted until `validate_known.py` passes. It rebuilds known
families from scratch and compares against the literature:

```
$ sage validate_known.py
1. Barreto-Naehrig (k=12, D=3, rho=1)
  [PASS] BN from Phi_12(6x^2) factorisation
  [PASS] BN-equivalent rho=1 found by dense sweep of Q(zeta_12)
2. Cyclotomic and KSS families (structured search w = zeta_ell * pi)
  [PASS] KSS16: rho = 5/4      [PASS] KSS18: rho = 4/3
  [PASS] KSS36: rho = 7/6      [PASS] SG54:  rho = 10/9
  [PASS] KSS8:  rho = 5/4      [PASS] BLS24: rho = 5/4
  [PASS] BLS12: rho = 3/2      [PASS] BLS48: rho = 9/8
3. Freeman-Scott-Teske Construction 6.7 at k=72 (D=2)
  [PASS] closed form 6.7 == direct Brezing-Weng with D=2
  [PASS] rho = 4/3 matches the (5k/12+2)/phi(k) formula
4. Gasnier-Guillevic subfield method, Example 3.4 (k=22, D=7)
  [PASS] R, T, Y, Q all match the published polynomials
13/13 checks passed
```

Several of these match *bit for bit*, not merely in ρ: KSS16's
`r = x^8+48x^4+625`, KSS36's `r = x^12+683x^6+117649`, Scott–Guillevic's
`r = 3^9u^18+3^5u^9+1`, BN's `r = 36x^4+36x^3+18x^2+6x+1`, and all four of
Gasnier–Guillevic's Example 3.4 polynomials.

Runtime: ~30 s. `sage validate_known.py quick` skips the k=48 case.

## The framework

All constructions here are instances of one recipe (Freeman–Scott–Teske §6).
Pick a number field `K` containing `zeta_k` and `sqrt(-D)`, pick a primitive
element `theta` of `K`, and set

```
R = minpoly(theta)
T = representative of zeta_k + 1               (mod R)
Y = representative of (zeta_k - 1)/sqrt(-D)    (mod R)
Q = (T^2 + D*Y^2)/4
rho = deg Q / deg R
```

Three of the five validity conditions are then automatic, which is what makes
a large sweep tractable: `T = zeta_k+1` forces `R | Phi_k(T-1)`, the CM
equation `4Q-T^2 = D*Y^2` holds identically, and `Q = zeta_k mod R` forces
`R | Q+1-T`. The only real filters are on `Q`: irreducible, positive leading
coefficient, and representing primes.

The choice of `theta` is the whole game, and the modules differ only in where
they look for it.

| module | `K` | `theta` ranges over |
|---|---|---|
| `fast.py` | `Q(zeta_ell)` | all primitive elements in a coefficient box (dense or sparse) — this is the KSS method |
| `structured.py` | `Q(zeta_ell)` | `zeta_ell * pi`, `pi` a small Gaussian/Eisenstein integer |
| `l3.py` | `Q[x]/(r)`, `r \| Phi_ell(u(x))` | the natural generator — the non-cyclotomic layer where BN lives |
| `subfield.py` | `F.C_k`, `F = Q(sqrt(-D))` | `alpha * zeta_k`, `alpha` in `O_F` — Gasnier–Guillevic |

### A note on `structured.py`

Every published KSS-type family is `theta = zeta_ell * pi` for a small
Gaussian or Eisenstein integer `pi`, which makes `theta^(n/2)` land in the
quadratic subfield and gives

```
R(x) = x^n - Tr(theta^(n/2)) x^(n/2) + N(pi)^(n/2)
```

explaining the constant terms of the published `R`: KSS16 `5^4`, KSS32 `13^8`,
KSS18 `7^3`, KSS36 `7^6`, SG54 `3^9`. This is *not* a new observation — it is
the KSS method, and Remark 3.5 of Gasnier–Guillevic states the same
`theta = (a*zeta_d + b)*zeta_k` shape. It is implemented here because it makes
the sweep cheap.

### A note on `subfield.py`

Gasnier–Guillevic's method, implemented in closed form so no large linear
algebra is needed. The conjugates of `theta = alpha*zeta_k` are
`alpha*zeta_k^j` and `conj(alpha)*zeta_k^j`, so

```
R(X) = [alpha^phi * Phi_k(X/alpha)] * [conj(alpha)^phi * Phi_k(X/conj(alpha))]
```

and `T`, `Y` follow from 2x2 systems in `F` via `T = P1(X^e)X + 1`,
`Y = P2(X^e)X - P3(X^e)`, where `e = k/2` for even `k` with
`sqrt(-D)` outside `Q(zeta_k)`.

## Files

**Engine**

- `fast.py` — dense/sparse sweep of primitive elements of `Q(zeta_ell)`
- `structured.py` — the `zeta_ell * pi` family
- `subfield.py` — Gasnier–Guillevic subfield method, closed form
- `l3.py` — non-cyclotomic `K` from factors of `Phi_ell(u(x))`
- `normalize.py` — reparametrisation `x -> L*u + b`, content normalisation,
  and seed verification. Essential and easy to get wrong: the standard KSS18
  family is `r(x) = (x^6+37x^3+343)/343` on `x = 14 mod 42`, so a raw minimal
  polynomial must have its content divided out before it can represent primes.
- `pplift.py` — integral residues modulo large prime powers, via `GF(p)`
  root-finding plus Hensel lifting. Needed because subfield-method
  denominators reach `2^43` and `2603047^2`, far past any residue scan.
- `verify.py` — gold-standard check: find integer seeds giving `q`, `r` prime
  with `r | q+1-t`, CM equation solvable, and `ord_r(q)` exactly `k`.
- `famcheck.py` — decides whether a *potential* family is an actual family
  (Bunyakovsky criterion). Also a CLI:

  ```
  $ sage famcheck.py 22 7 0 1
  k=22 D=7 alpha=0+1w -> FAMILY
    rho=6/5  deg R=20  deg Q=24
    Q integral and content-free on x = 7406*u + 4
    content(R) = 541696 (divide R by this)
  ```

**Drivers**

- `validate_known.py` — the validation suite above
- `run_sweep.py` — sweep `fast.py` over a range of `k`
- `run_struct.py` — sweep `structured.py` over a range of `ell`
- `run_subfield.py` — sweep `subfield.py` over `k`, `D`, `alpha`
- `hunt.py` — compares unit vs non-unit `pi` at fixed `k`

## Requirements

SageMath (developed on 10.8). No other dependencies.

```
sage validate_known.py
sage run_subfield.py 22,46 4 400
```

## References

- Freeman, Scott, Teske, *A taxonomy of pairing-friendly elliptic curves*,
  J. Cryptol. 23 (2010) 224–280.
- Kachisa, Schaefer, Scott, *Constructing Brezing-Weng pairing-friendly
  elliptic curves using elements in the cyclotomic field*, Pairing 2008.
- Scott, Guillevic, *A new family of pairing-friendly elliptic curves*,
  WAIFI 2018 (eprint 2018/193).
- Gasnier, Guillevic, *An algebraic point of view on the generation of
  pairing-friendly curves*, SIAM J. Appl. Algebra Geom. 9 (2025) 456–480.

## Licence

MIT.
