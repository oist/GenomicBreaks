"""Stage 1: every candidate on every ordered pair, on the alignment as it is
("fwd") and after exchanging target and query ("swp").  One npz per Branch."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
import numpy as np
import mk_common as mk, mk_engine as en

NAMES = en.candidate_names("raw") + en.candidate_names("coa")
OUT = os.path.join(mk.RES, "stage1")


def one(args):
    br, t, q = args
    t0 = time.time()
    b = mk.load_cache(br, t, q)
    b.pop("tnames"); b.pop("qnames")
    res = []
    for bb in (b, mk.swap(b)):
        f = en.features(bb, "raw")
        f.update(en.features(mk.coalesce(bb), "coa"))
        res.append(np.array([f[x] for x in NAMES], dtype=np.float64))
    return br, t, q, res[0], res[1], len(b["ts"]), time.time() - t0


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "names.txt"), "w") as fh:
        fh.write("\n".join(NAMES) + "\n")
    t0 = time.time()
    for br in mk.branches():
        fn = os.path.join(OUT, br + ".npz")
        if os.path.exists(fn):
            continue
        jobs = [(br, t, q) for t, q in mk.pairs_of(br)]
        fwd = np.full((len(jobs), len(NAMES)), np.nan)
        swp = np.full((len(jobs), len(NAMES)), np.nan)
        nblocks = np.zeros(len(jobs))
        index = {(t, q): i for i, (_, t, q) in enumerate(jobs)}
        with Pool(2) as pool:
            for _, t, q, a, b, n, dt in pool.imap_unordered(one, jobs, chunksize=1):
                i = index[(t, q)]
                fwd[i], swp[i], nblocks[i] = a, b, n
        np.savez(fn, fwd=fwd, swp=swp, nblocks=nblocks,
                 target=np.array([j[1] for j in jobs]), query=np.array([j[2] for j in jobs]))
        print(br, len(jobs), "pairs", int(nblocks.sum()), "blocks", round(time.time() - t0), "s", flush=True)
    print("finished", round(time.time() - t0), "s", flush=True)
