"""Tables read by the vignette of the R package (vignettes/MonkeyIndices.Rmd).

Written in the directory given as first argument:

  MonkeyIndices_values.tsv     one row per alignment: clade, genome pair, P
                               distance and the value of each selected function
  MonkeyIndices_functions.tsv  one row per selected function: its name in the
                               screen and its metrics
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import mk_common as mk

FUN = [("wukong_index", "raw.tailw.qstep"), ("abu_index", "raw.bowley.gdiff"), ("albert_index", "raw.kurt.gdiff"),
       ("marcel_index", "raw.sd.pdisp"), ("imo_index", "raw.xcor.tedge~qedge"), ("rafiki_index", "raw.rcor.qstep~jm_qedge"),
       ("hanuman_index", "raw.rcor.jm_tw~jd_qedge"), ("jack_index", "raw.rcor.rstep~dirc"),
       ("able_index", "raw.pcor.rstep~jd_qedge"), ("mizaru_index", "coa.rcor.tgap~jd_qedge"),
       ("george_index", "raw.pcor.tgap~jm_sdens"), ("baker_index", "raw.wcor.gapf~bgap")]
# Readable names of the reference statistics that come closest to a function
LABEL = {"index_synteny_query": "synteny index (query)", "index_synteny_target": "synteny index (target)",
         "index_strandRand_target": "strand randomisation index (target)",
         "index_strandRand_query": "strand randomisation index (query)",
         "bpg_dist_frac": "breakpoint-graph distance per alignment", "bpg_reuse": "breakpoint reuse",
         "guessed_query_length_Q1": "query sequence length (first quartile)",
         "translocated_width_target_Median": "translocated width (target, median)",
         "chain_width_target_L50": "chain width (target, L50)",
         "twin_inv_width_target_Min": "twin inversion width (target, minimum)"}

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
sig = lambda x, n: "" if not np.isfinite(x) else "%.*g" % (n, x)

fv = pd.read_csv(os.path.join(mk.RES, "final_values.tsv"), sep="\t")
key = [b + "|" + "|".join(sorted((t, q))) for b, t, q in zip(fv["branch"], fv["target"], fv["query"])]
val = pd.DataFrame({"clade": fv.branch.str.replace(r"^(MAM8|RY4)_", "", regex=True),
                    "pair": pd.factorize(np.array(key))[0] + 1,
                    "P": [sig(x, 5) for x in fv.P]})
for name, cand in FUN:
    val[name] = [sig(x, 5) for x in fv[cand]]
val.to_csv(os.path.join(out, "MonkeyIndices_values.tsv"), sep="\t", index=False)

met = pd.read_csv(os.path.join(mk.RES, "selected_final_metrics.tsv"), sep="\t", index_col=0)
rows = []
for name, cand in FUN:
    m = met.loc[cand]
    rows.append({"fun": name, "candidate": cand,
                 "median": sig(m.null_median_real, 6), "shuffled": sig(m.null_median_shuffled, 6),
                 "departs": sig(m.null_frac3, 4), "reciprocal": sig(m.recip_med, 4), "symmetry": sig(m.sym_med, 4),
                 "closest": m.ref_which, "closest_label": LABEL[m.ref_which], "closest_r": sig(m.ref_max, 4),
                 "others_r": sig(m.mutual_max, 4), "ridge_r": sig(m.tmr_ridge_r, 4), "trees_r": sig(m.tmr_trees_r, 4)})
pd.DataFrame(rows).to_csv(os.path.join(out, "MonkeyIndices_functions.tsv"), sep="\t", index=False)
print(len(val), "alignments,", val.pair.nunique(), "genome pairs,", len(rows), "functions")
