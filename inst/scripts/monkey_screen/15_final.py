"""Recompute the selected candidates with the final code on every pair, in
both orientations, and compare with the values used during the screen."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings("ignore")
from multiprocessing import Pool
import numpy as np, pandas as pd
import mk_common as mk, mk_engine as en, mk_screen as sc

WANT = pd.read_csv(os.path.join(mk.RES, sys.argv[1] if len(sys.argv) > 1 else "selected.tsv"),
                   sep="\t", index_col=0).index.tolist()
RAW = [w for w in WANT if w.startswith("raw.")]
COA = [w for w in WANT if w.startswith("coa.")]


def one(args):
    i, br, t, q = args
    b = mk.load_cache(br, t, q)
    b.pop("tnames"); b.pop("qnames")
    res = []
    for bb in (b, mk.swap(b)):
        f = {}
        if RAW:
            f.update(en.features(bb, "raw", want=RAW))
        if COA:
            f.update(en.features(mk.coalesce(bb), "coa", want=COA))
        res.append([f[w] for w in WANT])
    return i, res[0], res[1]


if __name__ == "__main__":
    t0 = time.time()
    design = sc.load_design()
    jobs = [(i, r.branch, r.target, r.query) for i, r in enumerate(design.itertuples())]
    fwd = np.full((len(jobs), len(WANT)), np.nan)
    swp = np.full((len(jobs), len(WANT)), np.nan)
    with Pool(2) as pool:
        for i, a, b in pool.imap_unordered(one, jobs, chunksize=4):
            fwd[i], swp[i] = a, b
    out = design[["branch", "target", "query", "P"]].copy()
    for k, w in enumerate(WANT):
        out[w] = fwd[:, k]
    for k, w in enumerate(WANT):
        out["swap:" + w] = swp[:, k]
    out.to_csv(os.path.join(mk.RES, "final_values.tsv"), sep="\t", index=False)
    names, f1, s1 = sc.load_stage1(design)
    col = [names.index(w) for w in WANT]
    for lab, new, old in (("fwd", fwd, f1[:, col]), ("swp", swp, s1[:, col])):
        both = np.isfinite(new) & np.isfinite(old)
        d = np.where(both, np.abs(new - old) / np.maximum(np.abs(old), 1e-9), 0)
        print(lab, "max relative difference with stage 1, per function:", np.round(d.max(axis=0), 8),
              "| defined/undefined changes:", int((np.isfinite(new) != np.isfinite(old)).sum()))
    print("done", round(time.time() - t0), "s")
