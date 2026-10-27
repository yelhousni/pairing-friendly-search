import os, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sage.all import euler_phi, QQ
from fast import sweep_fast

mode = sys.argv[1]

if mode == 'validate':
    for (k, ell, cr, mt, name, pub) in [(16,16,2,2,'KSS16','5/4'),
                                        (18,18,2,2,'KSS18','4/3'),
                                        (8, 24,2,2,'KSS8','5/4'),
                                        (36,36,2,2,'KSS36','7/6'),
                                        (54,54,1,2,'SG54','10/9'),
                                        (12,12,2,2,'BLS12','3/2'),
                                        (24,24,2,2,'BLS24','5/4')]:
        t0 = time.time()
        res = sweep_fast(k, ell, cr, mt, Dmax=12)
        if not res:
            print(f"{name}: none [{time.time()-t0:.1f}s]"); continue
        bD, bp = min(((D, v[0]) for D, v in res.items()), key=lambda p: p[1][0])
        print(f"{name:7s} k={k} ell={ell}: best rho={bp[0]} D={bD} "
              f"(published {pub}) [{time.time()-t0:.1f}s]")
        print(f"         r = {bp[1]}")
        print(f"         t = {bp[2]}")
    sys.exit()

KMIN, KMAX, NMAX = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
OUT = sys.argv[5]

def tier(n):
    if n <= 8:   return (3, 3, None)
    if n <= 12:  return (3, 3, None)
    if n <= 16:  return (2, 3, None)
    if n <= 20:  return (2, 2, None)
    return (2, 2, 16)

fout = open(OUT, 'a')
T0 = time.time()
for k in range(KMIN, KMAX + 1):
    for ell in range(k, 12 * k + 1, k):
        n = euler_phi(ell)
        if n < 2 or n > NMAX:
            continue
        cr, mt, mi = tier(n)
        t0 = time.time()
        try:
            res = sweep_fast(k, ell, cr, mt, Dmax=20, max_idx=mi)
        except Exception as e:
            print(f"k={k} ell={ell} ERR {e}", flush=True); continue
        if not res:
            continue
        for D, cands in res.items():
            for (rho, r, t, y, q, winfo) in cands[:2]:
                fout.write(json.dumps(dict(
                    k=int(k), ell=int(ell), n=int(n), D=int(D),
                    rho=str(rho), rho_f=float(rho),
                    w=str(winfo), r=str(r), t=str(t), y=str(y), q=str(q))) + "\n")
        fout.flush()
        br = min(v[0][0] for v in res.values())
        print(f"k={k:3d} ell={ell:4d} n={n:2d} best_rho={br} "
              f"[{time.time()-t0:.1f}s]", flush=True)
print(f"DONE {time.time()-T0:.1f}s", flush=True)
