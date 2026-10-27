"""
The decisive experiment.

For each field Q(zeta_ell) and each k | ell, compare:
  (A) the best rho over UNIT pi  (N(pi)=1)  -- this is the classical
      cyclotomic Brezing-Weng/BLS family, r = Phi_ell up to scaling;
  (B) the best rho over ALL pi               -- the KSS-type construction.

Published KSS-type gains occur exactly at ell in {16,18,32,36,54}. Any other
ell where (B) strictly beats (A), with a fully verified family, is a
construction outside the published set.
"""

import os, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sage.all import euler_phi, QQ, ZZ, divisors, gcd
from structured import search_ell
from normalize import verify_family, progressions, _denom


def norm_pi(D, a, b):
    a = ZZ(a); b = ZZ(b)
    if D % 4 == 3:
        return a * a + a * b + b * b * (1 + D) // 4
    return a * a + D * b * b


def run(LMIN, LMAX, NMAX, B, out_path, denom_cap=200000, tries=250,
        umax=1200):
    fout = open(out_path, 'a')
    for ell in range(LMIN, LMAX + 1):
        n = euler_phi(ell)
        if n < 4 or n > NMAX:
            continue
        ks = [d for d in divisors(ell) if d >= 5 and euler_phi(d) >= 4]
        if not ks:
            continue
        t0 = time.time()
        try:
            res = search_ell(ell, Dmax=20, B=B, k_list=ks, rho_cap=QQ(2))
        except Exception as e:
            print(f"ell={ell} ERR {e}", flush=True)
            continue
        if not res:
            continue
        for f in res:
            f['N'] = norm_pi(f['D'], f['a'], f['b'])
        groups = {}
        for f in res:
            groups.setdefault((f['k'], f['D']), []).append(f)

        for (k, D), cands in sorted(groups.items()):
            def best_verified(pool):
                pool = [f for f in pool if _denom(f['q']) <= denom_cap]
                pool.sort(key=lambda f: (f['rho'], f['N'],
                                         abs(f['a']) + abs(f['b'])))
                for f in pool[:tries]:
                    if not progressions(f['rI'], f['t'], f['q']):
                        continue
                    s = verify_family(f['rI'], f['t'], f['q'], D, k,
                                      umax=umax, want=1)
                    if s:
                        return f, s[0]
                return None, None

            fU, sU = best_verified([f for f in cands if f['N'] == 1])
            fA, sA = best_verified(cands)
            if fA is None:
                continue
            rhoU = fU['rho'] if fU else None
            # a real gain requires an actual verified baseline to beat
            gain = (rhoU is not None) and (fA['rho'] < rhoU)
            rec = dict(ell=int(ell), k=int(k), D=int(D), n=int(n),
                       rho_all=str(fA['rho']), rho_all_f=float(fA['rho']),
                       rho_unit=(str(rhoU) if rhoU else None),
                       N=int(fA['N']), a=int(fA['a']), b=int(fA['b']),
                       gain=bool(gain), r=str(fA['rI']), t=str(fA['t']),
                       q=str(fA['q']),
                       seed=dict(u0=int(sA['u0']), L=int(sA['L']),
                                 b=int(sA['b']), qbits=int(sA['qbits']),
                                 rbits=int(sA['rbits'])))
            fout.write(json.dumps(rec) + "\n")
            fout.flush()
            flag = "  <== NON-UNIT GAIN" if gain and fA['N'] > 1 else ""
            print(f"ell={ell:4d} k={k:4d} D={D:3d} rho_all={str(fA['rho']):7s}"
                  f" rho_unit={str(rhoU):7s} N(pi)={fA['N']:5d}{flag}",
                  flush=True)
        print(f"  -- ell={ell} done [{time.time()-t0:.1f}s]", flush=True)


if __name__ == '__main__':
    run(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]),
        int(sys.argv[4]), sys.argv[5])
