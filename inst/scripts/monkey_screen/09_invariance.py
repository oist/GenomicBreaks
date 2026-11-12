"""Representation test: the same alignments written with sequences reversed
at random and sequence codes permuted must give the same candidate values."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
import numpy as np, pandas as pd
import mk_common as mk, mk_engine as en, mk_null as nu, mk_screen as sc

NAMES = en.candidate_names("raw") + en.candidate_names("coa")


def one(args):
    i, br, t, q = args
    b = mk.load_cache(br, t, q)
    b.pop("tnames"); b.pop("qnames")
    rng = np.random.default_rng(1000 + i)
    out = []
    for bb in (b, nu.rerepresent(b, rng)):
        f = en.features(bb, "raw")
        f.update(en.features(mk.coalesce(bb), "coa"))
        out.append(np.array([f[x] for x in NAMES]))
    return i, out[0], out[1]


if __name__ == "__main__":
    design = sc.load_design()
    names, fwd, _ = sc.load_stage1(design)
    assert names == NAMES
    nb = design["aligned_width_target_N"].to_numpy()
    rng = np.random.default_rng(7)
    jobs = []
    for br in design.branch.unique():
        idx = np.flatnonzero((design.branch == br).to_numpy() & (nb <= 200000))
        for i in rng.choice(idx, min(20, len(idx)), replace=False):
            r = design.iloc[i]
            jobs.append((int(i), r.branch, r.target, r.query))
    worst = np.zeros(len(NAMES))
    nanflip = np.zeros(len(NAMES))
    with Pool(2) as pool:
        for i, a, v in pool.imap_unordered(one, jobs, chunksize=2):
            both = np.isfinite(a) & np.isfinite(v)
            d = np.where(both, np.abs(a - v) / np.maximum(np.maximum(np.abs(a), np.abs(v)), 1e-9), 0.0)
            # values that are numerically zero in both representations agree
            d[both & (np.abs(a) < 1e-9) & (np.abs(v) < 1e-9)] = 0.0
            worst = np.maximum(worst, d)
            nanflip += np.isfinite(a) != np.isfinite(v)
    out = pd.DataFrame({"inv_maxdiff": worst, "inv_nanflip": nanflip}, index=NAMES)
    out.index.name = "candidate"
    out.to_csv(os.path.join(mk.RES, "invariance.tsv"), sep="\t")
    print(len(jobs), "pairs tested")
    print("candidates with a relative difference > 1e-6:", int((worst > 1e-6).sum()),
          "; with a defined/undefined change:", int((nanflip > 0).sum()))
    print(out[(worst > 1e-6) | (nanflip > 0)].sort_values("inv_maxdiff", ascending=False).head(40).to_string())
