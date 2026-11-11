"""Linear transferable redundancy for all candidates."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import mk_common as mk, mk_screen as sc, mk_redundancy as rd

t0 = time.time()
design = sc.load_design()
names, fwd, swp = sc.load_stage1(design)
z = np.load(os.path.join(mk.RES, "ref_corr.npz"))
ref_cols = list(z["refs"])
E = design[ref_cols].to_numpy(float)
branch = design["branch"].to_numpy()
X = rd.within_rank(E, branch)
Y = rd.within_rank(fwd, branch)
print("ranks", round(time.time() - t0), flush=True)
r = rd.ridge(Y, X, branch)
out = pd.DataFrame({"tmr_ridge": r}, index=names)
out.index.name = "candidate"
out.to_csv(os.path.join(mk.RES, "redundancy_ridge.tsv"), sep="\t")
np.savez_compressed(os.path.join(mk.RES, "within_ranks.npz"), X=X, Y=Y)
print(out.describe().to_string(), round(time.time() - t0), "s")
