import os, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sage.all import (QQ, ZZ, euler_phi, is_squarefree, gcd, lcm,
                      PolynomialRing, is_prime)
from subfield import build, audit, rho_formula, sqrt_minus_D_in_Ck
from normalize import verify_family, progressions, _denom


def fst_best_fixedD(k):
    """Best FST fixed-discriminant rho we can state in closed form at this k.
    6.3: k = 2m, m odd  -> (k/2+2)/phi(k)   [D=1]
    6.4: k = 4m, m odd  -> (k/2+2)/phi(k)   [D=1]
    6.6: 18 not| k      -> (l/3 + c)/phi(l) [D=3], l = lcm(6,k),
                           c = 6 if k=4 mod 6 else 2
    6.7: 3 | k          -> (5k/12+2)/phi(k) [D=2] (k even)
    """
    phi = euler_phi(k)
    cands = []
    if k % 4 == 2 and (k // 2) % 2 == 1:
        cands.append((QQ(k // 2 + 2) / QQ(phi), '6.3'))
    if k % 4 == 0 and (k // 4) % 2 == 1:
        cands.append((QQ(k // 2 + 2) / QQ(phi), '6.4'))
    if k % 18 != 0:
        l = lcm(6, k)
        c = 6 if k % 6 == 4 else 2
        cands.append((QQ(l // 3 + c) / QQ(euler_phi(l)), '6.6'))
    if k % 3 == 0 and k % 2 == 0:
        cands.append((QQ(5 * k // 12 + 2) / QQ(phi), '6.7'))
    if not cands:
        return (None, None)
    return min(cands)


def search_k(k, Dmax=23, B=5, umax=600, denom_cap=10 ** 14, want_all=False):
    """Find an ACTUAL subfield-method family at embedding degree k.
    Candidates are tried smallest-first (small D, small alpha), since those
    give the smallest denominators and hence the best chance that Q is
    integral somewhere."""
    from normalize import progressions_crt
    Ds = [D for D in range(1, Dmax + 1)
          if is_squarefree(D) and not sqrt_minus_D_in_Ck(k, D)]
    cands = []
    for D in Ds:
        for a in range(-B, B + 1):
            for b in range(-B, B + 1):
                if b == 0 or gcd(a, b) != 1:
                    continue
                cands.append((abs(a) + abs(b), D, a, b))
    cands.sort()
    found = []
    for (_, D, a, b) in cands:
        try:
            f = build(k, D, a, b)
        except Exception:
            continue
        if f is None or audit(f):
            continue
        if _denom(f['Q']) > denom_cap:
            continue
        if not progressions_crt(f['R'], f['T'], f['Q']):
            continue
        s = verify_family(f['R'], f['T'], f['Q'], D, k, umax=umax, want=1,
                          use_crt=True)
        if s:
            f['seed'] = s[0]
            if not want_all:
                return f
            found.append(f)
    return found if want_all else None


if __name__ == '__main__':
    ks = [int(v) for v in sys.argv[1].split(',')]
    B = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    umax = int(sys.argv[3]) if len(sys.argv) > 3 else 600
    for k in ks:
        phi = euler_phi(k)
        rs = rho_formula(k, 0)
        fb, which = fst_best_fixedD(k)
        gain = (fb is not None and rs < fb)
        print(f"\n### k={k} phi={phi}  subfield rho={rs} ({float(rs):.4f})  "
              f"FST best={fb} via {which}  {'IMPROVES' if gain else ''}",
              flush=True)
        t0 = time.time()
        f = search_k(k, B=B, umax=umax)
        if f is None:
            print(f"    no verified family found [{time.time()-t0:.0f}s]", flush=True)
            continue
        s = f['seed']
        print(f"    FOUND D={f['D']} alpha={f['a']}{f['b']:+d}w  rho={f['rho']} "
              f"[{time.time()-t0:.0f}s]", flush=True)
        print(f"    deg R={f['R'].degree()} deg Q={f['Q'].degree()} e={f['e']}", flush=True)
        print(f"    seed: u0={s['u0']} (x={s['L']}u+{s['b']}) "
              f"q={s['qbits']}b r={s['rbits']}b rho_emp={s['rho']:.4f}", flush=True)
        json.dump(dict(k=int(k), D=int(f['D']), a=int(f['a']), b=int(f['b']),
                       rho=str(f['rho']), R=str(f['R']), T=str(f['T']),
                       Q=str(f['Q']), Y=str(f['Y']),
                       seed_q_bits=int(s['qbits']), seed_r_bits=int(s['rbits'])),
                  open(f'sub_k{k}.json', 'w'))
