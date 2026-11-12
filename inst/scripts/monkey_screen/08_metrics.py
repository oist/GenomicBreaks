"""Screening metrics of all candidates that do not need extra computation:
distribution, reciprocal agreement, ancestor-node agreement, and rank
correlation with the reference (existing) statistics."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import mk_common as mk, mk_screen as sc

t0 = time.time()
design = sc.load_design()
names, fwd, swp = sc.load_stage1(design)
branch = design["branch"].to_numpy()
meta = ["branch", "target", "query", "mate", "Pm", "node", "upair"]
ref_cols = [c for c in design.columns if c not in meta and c != "P_matrix"
            and np.issubdtype(design[c].dtype, np.number)]
E = design[ref_cols].to_numpy(float)
# drop reference columns that are constant or almost never defined
keep = [j for j in range(E.shape[1]) if np.isfinite(E[:, j]).mean() >= 0.2 and np.nanstd(E[:, j]) > 0]
ref_cols = [ref_cols[j] for j in keep]
E = E[:, keep]
print(len(names), "candidates;", len(ref_cols), "reference statistics;", len(design), "pairs", flush=True)

M = {}
M.update(sc.distribution_metrics(fwd, branch))
print("distribution", round(time.time() - t0), flush=True)
M.update(sc.reciprocal_metrics(fwd, swp, design))
print("reciprocal", round(time.time() - t0), flush=True)
M.update(sc.node_metrics(fwd, design))
print("nodes", round(time.time() - t0), flush=True)
pooled, pooled_p, within, within_p = sc.corr_blocks(fwd, E, design)
print("references", round(time.time() - t0), flush=True)
np.savez(os.path.join(mk.RES, "ref_corr.npz"), pooled=pooled, pooled_p=pooled_p, within=within,
         within_p=within_p, names=np.array(names), refs=np.array(ref_cols))
ip = ref_cols.index("P")
M["rho_P_pool"] = pooled[:, ip]
M["rho_P_within"] = within[:, ip]
stack = np.abs(np.stack([pooled, pooled_p, within, within_p]))       # 4 x cand x ref
stack = np.where(np.isfinite(stack), stack, 0.0)
for k, lab in enumerate(("pool", "poolP", "within", "withinP")):
    M["ref_" + lab] = stack[k].max(axis=1)
best = stack.max(axis=0)
M["ref_max"] = best.max(axis=1)
M["ref_which"] = np.array(ref_cols)[best.argmax(axis=1)]
out = pd.DataFrame(M, index=names)
out.index.name = "candidate"
out.to_csv(os.path.join(mk.RES, "metrics.tsv"), sep="\t")
print(out.describe().T.to_string())
print("done", round(time.time() - t0), flush=True)
