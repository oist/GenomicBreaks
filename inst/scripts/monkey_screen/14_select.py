"""Selection of the functions (10 by default; the number is the first argument).

Fixed filters:
  defined        finite in >= 95 % of pairs, and in >= 90 % of the pairs of
                 every Branch but one
  not uniform    >= 95 % distinct values, most frequent value <= 2 %
  invariant      same value when sequences are reversed or renamed
  reciprocal     rank correlation between f(T___Q) and f(swap(Q___T)):
                 median over Branches >= 0.9, first quartile >= 0.8
  nodes          pairs sharing an ancestor node are more alike than chance
                 (z >= 5), also after removing the rank of P (z >= 3)

Two requirements are relaxed step by step until enough functions remain.
  tau   common ceiling, from 0.2 upwards by 0.1, on
          - |rank correlation| with every reference statistic (pooled,
            within Branch, each also partial on P);
          - multivariate transferable redundancy sqrt(R2), ridge then trees;
          - |rank correlation| with the function's own shuffled value;
          - |rank correlation| between selected functions.
  null  share of the genome pairs in which the function departs from 5
        shuffles by |z| >= 3: three quarters, then one half.
tau is raised only when the null requirement of one half is not enough:
correlation with existing statistics is the explicit request, the share of
pairs is a setting of this screen.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import mk_common as mk, mk_screen as sc, mk_redundancy as rd

R = mk.RES
design = sc.load_design()
names, fwd, swp = sc.load_stage1(design)
col = {n: i for i, n in enumerate(names)}
m = pd.read_csv(R + "/metrics_all.tsv", sep="\t", index_col=0).join(
    pd.read_csv(R + "/null_metrics.tsv", sep="\t", index_col=0))
m["null_r"] = m[["null_r_within", "null_r_pool"]].max(axis=1)
fixed = ((m.fin >= 0.95) & (m.fin_2nd >= 0.9) & (m.distinct >= 0.95) & (m["mode"] <= 0.02)
         & (m.inv_maxdiff <= 1e-4) & (m.inv_nanflip == 0)
         & (m.recip_med >= 0.9) & (m.recip_q25 >= 0.8)
         & (m.node_z >= 5) & (m.nodeR_z >= 3))
steps = [("defined", (m.fin >= 0.95) & (m.fin_2nd >= 0.9)),
         ("not uniform", (m.distinct >= 0.95) & (m["mode"] <= 0.02)),
         ("invariant", (m.inv_maxdiff <= 1e-4) & (m.inv_nanflip == 0)),
         ("reciprocal", (m.recip_med >= 0.9) & (m.recip_q25 >= 0.8)),
         ("nodes", (m.node_z >= 5) & (m.nodeR_z >= 3)),
         ("null, half of the pairs", m.null_frac3 >= 0.5),
         ("null, three quarters", m.null_frac3 >= 0.75)]
acc = pd.Series(True, index=m.index)
funnel = [("all candidates", len(m))]
for lab, f in steps:
    acc &= f
    funnel.append((lab, int(acc.sum())))
print("\n".join("%-16s %5d" % x for x in funnel))

pool = m[fixed & (m.null_frac3 >= 0.5) & (m.ref_max < 0.8)].index.tolist()
print(len(pool), "candidates pass the fixed filters with novelty < 0.8")
V = fwd[:, [col[n] for n in pool]]
mut = np.max(np.abs(np.nan_to_num(np.stack(sc.corr_blocks(V, V, design)))), axis=0)
mut = pd.DataFrame(np.maximum(mut, mut.T), index=pool, columns=pool)

zz = np.load(R + "/within_ranks.npz")
X, Y = zz["X"], zz["Y"]
branch = design["branch"].to_numpy()
tree_file = R + "/tmr_trees.json"
tree = json.load(open(tree_file)) if os.path.exists(tree_file) else {}


def tmr_trees(name):
    if name not in tree:
        tree[name] = rd.trees(Y[:, col[name]], X, branch)
        json.dump(tree, open(tree_file, "w"))
    return float(np.sqrt(tree[name]))


LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 10      # number of functions wanted


def greedy(tau, nullmin, limit=LIMIT):
    c = m.loc[pool]
    c = c[(c.null_frac3 >= nullmin) & (c.ref_max < tau) & (c.tmr_r < tau) & (c.null_r < tau)].copy()
    c["score"] = (c.recip_q25.rank(pct=True) + c.nodeR_eta.rank(pct=True)) / 2
    c = c.sort_values("score", ascending=False)
    picked, log = [], []
    for n in c.index:
        if picked and mut.loc[n, picked].max() >= tau:
            log.append((n, "correlated with " + mut.loc[n, picked].idxmax()))
            continue
        t = tmr_trees(n)
        if t >= tau:
            log.append((n, "predictable from references (trees) %.2f" % t))
            continue
        picked.append(n)
        if len(picked) == limit:
            break
    return picked, c, log


result = None
grid = []
for tau in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8):
    for nullmin in (0.75, 0.5):
        picked, c, log = greedy(tau, nullmin)
        grid.append((tau, nullmin, len(c), len(picked)))
        print("tau = %.1f, null >= %.2f: %d candidates before the mutual and tree filters, %d selected"
              % (tau, nullmin, len(c), len(picked)), flush=True)
        if len(picked) == LIMIT and result is None:
            result = (tau, nullmin, picked, c, log)
    if result is not None and tau >= result[0] + 0.1 - 1e-9:
        break

tau, nullmin, picked, c, log = result
sel = c.loc[picked].copy()
sel["tmr_trees_r"] = [tmr_trees(n) for n in picked]
sel["mutual_max"] = [mut.loc[n, [p for p in picked if p != n]].max() for n in picked]
sel.to_csv(R + "/selected.tsv", sep="\t")
c.to_csv(R + "/passing_at_tau.tsv", sep="\t")
pd.DataFrame(log, columns=["candidate", "reason"]).to_csv(R + "/skipped_before_last.tsv", sep="\t", index=False)
mut.loc[picked, picked].to_csv(R + "/selected_mutual.tsv", sep="\t")
json.dump({"tau": tau, "nullmin": nullmin, "funnel": funnel, "pool": len(pool), "passing_at_tau": len(c),
           "grid": grid},
          open(R + "/selection_summary.json", "w"))
pd.set_option("display.width", 250)
print(sel[["score", "recip_med", "recip_q25", "sym_med", "node_eta", "nodeR_eta", "genome_r2", "null_frac3",
           "null_medz", "null_r", "rho_P_within", "ref_max", "ref_which", "tmr_r", "tmr_trees_r", "mutual_max"]].round(3).to_string())
