"""Candidates for which the shuffling null is computed."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import mk_common as mk, mk_screen as sc, mk_engine as en

R = mk.RES
design = sc.load_design()
names, fwd, swp = sc.load_stage1(design)
m = (pd.read_csv(R + "/metrics.tsv", sep="\t", index_col=0)
     .join(pd.read_csv(R + "/invariance.tsv", sep="\t", index_col=0))
     .join(pd.read_csv(R + "/redundancy_ridge.tsv", sep="\t", index_col=0)))
m["tmr_r"] = np.sqrt(m.tmr_ridge)
branch = design["branch"].to_numpy()
fb = np.array([np.isfinite(fwd[branch == b]).mean(axis=0) for b in pd.unique(branch)])
m["fin_2nd"] = np.sort(fb, axis=0)[1]            # second lowest Branch
m["free"] = [en.arrangement_free(n) for n in m.index]
m["level"] = [n.split(".")[0] for n in m.index]
m["reducer"] = [n.split(".")[1] for n in m.index]
m.to_csv(R + "/metrics_all.tsv", sep="\t")
base = ((m.fin >= 0.95) & (m.fin_2nd >= 0.9) & (m.distinct >= 0.95) & (m["mode"] <= 0.02)
        & (m.inv_maxdiff <= 1e-4) & (m.inv_nanflip == 0))
pre = base & (m.recip_med >= 0.8) & (m.recip_q25 >= 0.7) & (m.node_z >= 5) & (m.nodeR_z >= 3)
print("defined, not uniform, invariant:", int(base.sum()), "| + reciprocal and nodes (relaxed 0.8/0.7):", int(pre.sum()))
for t in (0.5, 0.6, 0.7, 0.8):
    s = pre & (m.ref_max < t) & (m.tmr_r < t)
    print("  novelty <", t, ":", int(s.sum()), "of which arrangement-free:", int((s & m.free).sum()),
          "| raw", int((s & ~m.free & (m.level == "raw")).sum()), "coa", int((s & ~m.free & (m.level == "coa")).sum()))
short = m[pre & (m.ref_max < 0.7) & (m.tmr_r < 0.7) & ~m.free]
with open(R + "/shortlist.txt", "w") as fh:
    fh.write("\n".join(short.index) + "\n")
print(len(short), "candidates in the shortlist")
print(short.groupby(["level", "reducer"]).size().unstack(0).fillna(0).astype(int).T.to_string())
