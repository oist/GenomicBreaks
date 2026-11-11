"""All screening metrics of the selected functions, recomputed from the
values given by the final code."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import mk_common as mk, mk_screen as sc, mk_redundancy as rd

R = mk.RES
design = sc.load_design()
fv = pd.read_csv(R + "/final_values.tsv", sep="\t")
assert (fv.target.values == design.target.values).all() and (fv["query"].values == design["query"].values).all()
cand = [c for c in fv.columns if c.startswith("raw.") or c.startswith("coa.")]
fwd = fv[cand].to_numpy(float)
swp = fv[["swap:" + c for c in cand]].to_numpy(float)
branch = design["branch"].to_numpy()
M = {}
M.update(sc.distribution_metrics(fwd, branch))
fb = np.array([np.isfinite(fwd[branch == b]).mean(axis=0) for b in pd.unique(branch)])
M["fin_2nd"] = np.sort(fb, axis=0)[1]
M.update(sc.reciprocal_metrics(fwd, swp, design))
M.update(sc.node_metrics(fwd, design))
z = np.load(R + "/ref_corr.npz")
ref_cols = list(z["refs"])
E = design[ref_cols].to_numpy(float)
blocks = sc.corr_blocks(fwd, E, design)
stack = np.abs(np.nan_to_num(np.stack(blocks)))
best = stack.max(axis=0)
M["ref_max"] = best.max(axis=1)
M["ref_which"] = np.array(ref_cols)[best.argmax(axis=1)]
M["ref_variant"] = np.array(["pooled", "pooled, partial on P", "within Branch", "within Branch, partial on P"])[
    stack[:, np.arange(len(cand)), best.argmax(axis=1)].argmax(axis=0)]
ip = ref_cols.index("P")
M["rho_P_pool"] = blocks[0][:, ip]
M["rho_P_within"] = blocks[2][:, ip]
X = rd.within_rank(E, branch)
Y = rd.within_rank(fwd, branch)
M["tmr_ridge_r"] = np.sqrt(rd.ridge(Y, X, branch))
M["tmr_trees_r"] = np.sqrt([rd.trees(Y[:, j], X, branch) for j in range(len(cand))])
# null
names = [l.strip() for l in open(R + "/stage1/names.txt")]
col = [names.index(c) for c in cand]
key = {(r.branch, r.target, r.query): i for i, r in enumerate(design.itertuples())}
rows, nulls = [], []
for br in design.branch.unique():
    zz = np.load(os.path.join(R, "null", br + ".npz"))
    rows += [key[(br, t, q)] for t, q in zip(zz["target"], zz["query"])]
    nulls.append(zz["null"][:, :, col])
rows = np.array(rows)
null = np.concatenate(nulls)
real = fwd[rows]
mu, sd = np.nanmean(null, axis=1), np.nanstd(null, axis=1, ddof=1)
ok = np.isfinite(real) & np.isfinite(mu) & np.isfinite(sd)
zs_ = np.where(ok, (real - mu) / np.where(sd > 0, sd, np.nan), np.nan)
M["null_frac3"] = np.nanmean(np.abs(zs_) >= 3, axis=0)
M["null_medz"] = np.nanmedian(np.abs(zs_), axis=0)
M["null_median_real"] = np.nanmedian(real, axis=0)
M["null_median_shuffled"] = np.nanmedian(mu, axis=0)
rr = [sc.spearman_cols(real[branch[rows] == b], mu[branch[rows] == b]) for b in pd.unique(branch[rows])]
M["null_r"] = np.maximum(np.abs(np.nanmedian(np.array(rr), axis=0)), np.abs(sc.spearman_cols(real, mu)))
mut = np.max(np.abs(np.nan_to_num(np.stack(sc.corr_blocks(fwd, fwd, design)))), axis=0)
mut = np.maximum(mut, mut.T)
np.fill_diagonal(mut, 0)
M["mutual_max"] = mut.max(axis=1)
for q in (0.05, 0.25, 0.5, 0.75, 0.95):
    M["q%02d" % int(q * 100)] = np.nanquantile(fwd, q, axis=0)
out = pd.DataFrame(M, index=cand)
out.index.name = "candidate"
out.to_csv(R + "/selected_final_metrics.tsv", sep="\t")
pd.DataFrame(mut, index=cand, columns=cand).to_csv(R + "/selected_final_mutual.tsv", sep="\t")
# every existing reference above 0.4 for each function, for the report
top = []
for j, c in enumerate(cand):
    o = np.argsort(-best[j])[:4]
    top.append((c, "; ".join("%s %.2f" % (ref_cols[k], best[j, k]) for k in o)))
pd.DataFrame(top, columns=["candidate", "closest references"]).to_csv(R + "/selected_closest_refs.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(out[["fin", "fin_2nd", "distinct", "recip_med", "recip_q25", "sym_med", "node_z", "nodeR_z", "node_eta", "nodeR_eta", "genome_r2"]].round(3).to_string())
print(out[["null_frac3", "null_medz", "null_median_real", "null_median_shuffled", "null_r", "rho_P_within", "ref_max", "ref_which", "tmr_ridge_r", "tmr_trees_r", "mutual_max"]].round(3).to_string())
print(pd.DataFrame(top).to_string())
