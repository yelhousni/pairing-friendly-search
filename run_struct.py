import os, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sage.all import euler_phi, QQ, ZZ, divisors, lcm
from structured import search_ell, rigorous_check

LMIN, LMAX, NMAX, B, OUT = (int(sys.argv[1]), int(sys.argv[2]),
                            int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
LCAP = 10 ** 7          # skip candidates whose q has an enormous denominator
MAX_TRIES = 80          # validation attempts per (k,D) rho-class


def denom(f):
    L = ZZ(1)
    for c in f.coefficients():
        L = lcm(L, QQ(c).denominator())
    return L


fout = open(OUT, 'a')
best = {}
T0 = time.time()
for ell in range(LMIN, LMAX + 1):
    n = euler_phi(ell)
    if n < 2 or n > NMAX:
        continue
    ks = [d for d in divisors(ell) if d >= 5 and euler_phi(d) >= 4]
    if not ks:
        continue
    t0 = time.time()
    try:
        res = search_ell(ell, Dmax=40, B=B, k_list=ks, rho_cap=QQ(2))
    except Exception as e:
        print(f"ell={ell} ERR {e}", flush=True)
        continue
    if not res:
        continue
    groups = {}
    for f in res:
        groups.setdefault((f['k'], f['D']), []).append(f)
    nvalid = 0
    for (k, D), cands in groups.items():
        # cheapest first: small rho, then small denominator, then small pi
        cands.sort(key=lambda f: (f['rho'], abs(f['a']) + abs(f['b'])))
        tries = 0
        for f in cands:
            if tries >= MAX_TRIES:
                break
            if denom(f['q']) > LCAP:
                continue
            tries += 1
            if rigorous_check(f['rI'], f['t'], f['y'], f['q'], k, D):
                continue
            nvalid += 1
            if k not in best or f['rho'] < best[k]['rho']:
                best[k] = f
            fout.write(json.dumps(dict(
                k=int(k), ell=int(ell), D=int(D), n=int(n),
                a=int(f['a']), b=int(f['b']),
                rho=str(f['rho']), rho_f=float(f['rho']),
                r=str(f['rI']), t=str(f['t']), q=str(f['q']))) + "\n")
            break          # best valid one for this (k,D) at this ell
    fout.flush()
    print(f"ell={ell:4d} n={n:2d} cands={len(res):6d} groups={len(groups):3d} "
          f"valid={nvalid:3d} [{time.time()-t0:.1f}s]", flush=True)

print("\n=== BEST rho per k ===", flush=True)
for k in sorted(best):
    f = best[k]
    print(f"k={k:4d} rho={str(f['rho']):8s} D={f['D']:3d} ell={f['ell']:4d} "
          f"pi={f['a']}+{f['b']}w  r={f['rI']}", flush=True)
print(f"TOTAL {time.time()-T0:.1f}s", flush=True)
