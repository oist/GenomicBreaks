"""Departure of every candidate from its shuffled value."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import mk_common as mk, mk_screen as sc

design = sc.load_design()
names, fwd, swp = sc.load_stage1(design)
key = {(r.branch, r.target, r.query): i for i, r in enumerate(design.itertuples())}
rows, nulls = [], []
for br in design.branch.unique():
    fn = os.path.join(mk.RES, "null", br + ".npz")
    if not os.path.exists(fn):
        print("WARNING: no null for", br)
        continue
    z = np.load(fn)
    rows += [key[(br, t, q)] for t, q in zip(z["target"], z["query"])]
    nulls.append(z["null"])
rows = np.array(rows)
null = np.concatenate(nulls)                      # pairs x K x candidates
real = fwd[rows]
branch = design["branch"].to_numpy()[rows]
mu = np.nanmean(null, axis=1)
sd = np.nanstd(null, axis=1, ddof=1)
ok = np.isfinite(real) & np.isfinite(mu) & np.isfinite(sd)
with np.errstate(all="ignore"):
    diff = real - mu
    # differences at the level of floating-point noise are no difference
    same = np.abs(diff) <= 1e-8 * np.maximum(np.maximum(np.abs(real), np.abs(mu)), 1e-12)
    zz = np.where(same, 0.0, np.where(sd > 0, diff / sd, np.inf))
    zz = np.where(ok, zz, np.nan)
    out = pd.DataFrame(index=names)
    out.index.name = "candidate"
    out["null_n"] = ok.sum(axis=0)
    out["null_frac3"] = np.nanmean(np.abs(zz) >= 3, axis=0)
    out.loc[out.null_n == 0, "null_frac3"] = np.nan
    out["null_medz"] = np.nanmedian(np.abs(zz), axis=0)
    q1, q3 = np.nanquantile(real, [0.25, 0.75], axis=0)
    out["null_shift"] = np.nanmedian(np.abs(real - mu), axis=0) / (q3 - q1)
    # is the real value predictable from the shuffled value?
    rr = []
    for b in pd.unique(branch):
        m = branch == b
        rr.append(sc.spearman_cols(np.where(ok[m], real[m], np.nan), np.where(ok[m], mu[m], np.nan)))
    out["null_r_within"] = np.abs(np.nanmedian(np.array(rr), axis=0))
    out["null_r_pool"] = np.abs(sc.spearman_cols(np.where(ok, real, np.nan), np.where(ok, mu, np.nan)))
out = out.fillna({"null_r_within": 0.0, "null_r_pool": 0.0})
out.to_csv(os.path.join(mk.RES, "null_metrics.tsv"), sep="\t")
print(len(rows), "pairs;", out.describe().T.to_string())
