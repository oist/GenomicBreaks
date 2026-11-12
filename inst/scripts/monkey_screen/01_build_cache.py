"""Parse every blasttabplus file once and store the raw blocks as npz."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
import numpy as np
import mk_common as mk


def one(args):
    branch, t, q = args
    out = mk.cache_file(branch, t, q)
    if os.path.exists(out):
        return (branch, t, q, -1)
    b = mk.read_blasttabplus(mk.path(branch, "blasttabplus", t, q))
    mk.save_cache(branch, t, q, b)
    return (branch, t, q, len(b["ts"]))


if __name__ == "__main__":
    jobs = [(br, t, q) for br in mk.branches() for (t, q) in mk.pairs_of(br)]
    print(len(jobs), "pairs", flush=True)
    t0 = time.time()
    n = 0
    with Pool(2) as pool:
        for i, r in enumerate(pool.imap_unordered(one, jobs, chunksize=4), 1):
            n += max(r[3], 0)
            if i % 100 == 0:
                print(i, "done", n, "rows", round(time.time() - t0), "s", flush=True)
    print("finished", n, "rows", round(time.time() - t0), "s", flush=True)
