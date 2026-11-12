"""Small multiples: each selected function against divergence (P)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mk_common as mk

SURFACE, INK, QUIET, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
GROUPS = [("Prokaryotes (5 clades)", "#2a78d6", ["Buchnera", "Rickettsiales", "Neisseriales", "Deinococcota", "Nanobdellati"]),
          ("Fungi and Discoba (3 clades)", "#eb6834", ["Saccharomycotina", "Pezizomycotina", "Discoba"]),
          ("Animals (8 clades)", "#1baf7a", ["Apoidea", "Gastropoda", "Bivalvia", "Otomorpha", "Eupercaria", "Ovalentaria",
                                             "Other_Neoteleostei", "Primates"])]
FUN = [("wukong_index()", "raw.tailw.qstep"), ("abu_index()", "raw.bowley.gdiff"), ("albert_index()", "raw.kurt.gdiff"),
       ("marcel_index()", "raw.sd.pdisp"), ("imo_index()", "raw.xcor.tedge~qedge"), ("rafiki_index()", "raw.rcor.qstep~jm_qedge"),
       ("hanuman_index()", "raw.rcor.jm_tw~jd_qedge"), ("jack_index()", "raw.rcor.rstep~dirc"),
       ("able_index()", "raw.pcor.rstep~jd_qedge"), ("mizaru_index()", "coa.rcor.tgap~jd_qedge"),
       ("george_index()", "raw.pcor.tgap~jm_sdens"), ("baker_index()", "raw.wcor.gapf~bgap")]

fv = pd.read_csv(os.path.join(mk.RES, "final_values.tsv"), sep="\t")
met = pd.read_csv(os.path.join(mk.RES, "selected_final_metrics.tsv"), sep="\t", index_col=0)
fv["clade"] = fv.branch.str.replace("MAM8_", "").str.replace("RY4_", "")
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": GRID, "axes.linewidth": 0.6,
                     "xtick.color": QUIET, "ytick.color": QUIET, "text.color": INK, "axes.labelcolor": QUIET,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE})
fig, axes = plt.subplots(3, 4, figsize=(12.5, 9.6), sharex=True)
rng = np.random.default_rng(0)
order = rng.permutation(len(fv))                       # no group systematically on top
for ax, (name, cand) in zip(axes.ravel(), FUN):
    y = fv[cand].to_numpy()
    lo, hi = np.nanquantile(y, [0.01, 0.99])
    pad = 0.06 * (hi - lo)
    col = np.empty(len(fv), dtype=object)
    for _, c, clades in GROUPS:
        col[fv.clade.isin(clades).to_numpy()] = c
    ax.scatter(fv.P.to_numpy()[order], y[order], s=7, c=list(col[order]), linewidths=0.25, edgecolors=SURFACE, alpha=0.85, zorder=3)
    s = met.loc[cand, "null_median_shuffled"]
    if lo - pad <= s <= hi + pad:
        ax.axhline(s, color=QUIET, lw=0.8, ls=(0, (4, 3)), zorder=2)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlim(-0.01, 0.56)
    ax.set_title(name, loc="left", fontsize=9.5, fontweight="bold", color=INK, pad=5)
    ax.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.tick_params(length=0, labelsize=7.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
for ax in axes[-1]:
    ax.set_xlabel("P distance", fontsize=8)
fig.suptitle("The twelve monkey indices against divergence, one dot per alignment (1,910 alignments, 16 clades)",
             x=0.04, y=0.988, ha="left", fontsize=11.5, fontweight="bold", color=INK)
fig.text(0.04, 0.962, "Vertical axes cover the 1st to 99th percentiles of each index.  Dashed line: median of the index after shuffling, when it falls within the axis.",
         fontsize=8, color=QUIET, ha="left")
handles = [plt.Line2D([], [], marker="o", ls="", ms=5.5, mfc=c, mec=SURFACE, mew=0.4, label=lab) for lab, c, _ in GROUPS]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.033, 0.957), ncol=3, frameon=False, fontsize=8.5, handletextpad=0.3, columnspacing=1.6)
fig.subplots_adjust(left=0.05, right=0.99, top=0.885, bottom=0.05, wspace=0.26, hspace=0.2)
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(mk.ROOT, "deliver", "monkey_indices_vs_P")
os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
fig.savefig(out + ".pdf")
fig.savefig(out + ".png", dpi=140)
print("written", out)
