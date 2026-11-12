"""Shuffling null for all candidates: K random arrangements of each
alignment (one order per genome pair).  Raw-level candidates: the raw
blocks are shuffled.  Coalesced-level candidates: the coalesced blocks are
shuffled, and are not coalesced again."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
import numpy as np
import mk_common as mk, mk_engine as en, mk_null as nu

K = 5
NAMES = en.candidate_names("raw") + en.candidate_names("coa")
OUT = os.path.join(mk.RES, "null")


def one(args):
    i, br, t, q = args
    b = mk.load_cache(br, t, q)
    b.pop("tnames"); b.pop("qnames")
    c = mk.coalesce(b)
    rng = np.random.default_rng(50000 + i)
    res = np.full((K, len(NAMES)), np.nan)
    for k in range(K):
        f = en.features(nu.shuffle_both(b, rng), "raw")
        f.update(en.features(nu.shuffle_both(c, rng), "coa"))
        res[k] = [f[w] for w in NAMES]
    return t, q, res


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    n = 0
    for br in mk.branches():
        fn = os.path.join(OUT, br + ".npz")
        if os.path.exists(fn):
            continue
        jobs = [(n + k, br, t, q) for k, (t, q) in enumerate(p for p in mk.pairs_of(br) if p[0] < p[1])]
        n += 1000
        null = np.full((len(jobs), K, len(NAMES)), np.nan)
        index = {(j[2], j[3]): k for k, j in enumerate(jobs)}
        with Pool(2) as pool:
            for t, q, r in pool.imap_unordered(one, jobs, chunksize=1):
                null[index[(t, q)]] = r
        np.savez(fn, null=null, target=np.array([j[2] for j in jobs]), query=np.array([j[3] for j in jobs]))
        print(br, len(jobs), "pairs", round(time.time() - t0), "s", flush=True)
    print("finished", round(time.time() - t0), "s", flush=True)
