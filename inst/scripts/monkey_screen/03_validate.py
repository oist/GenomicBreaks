"""Compare the Python loader / coalescing with the YAML statistics."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import mk_common as mk

ex = pd.read_csv(mk.RES + "/existing_yaml.tsv", sep="\t").set_index(["branch", "target", "query"])
rng = np.random.default_rng(0)
idx = rng.choice(len(ex), 300, replace=False)
keys = ["aligned_width_target_N", "aligned_width_target_Total", "aligned_width_query_Total",
        "chain_width_target_N", "chain_width_target_Total", "chain_width_query_Total",
        "guessed_target_length_N", "guessed_query_length_N"]
worst = {k: 0.0 for k in keys}
bad = []
for i in idx:
    br, t, q = ex.index[i]
    b = mk.load_cache(br, t, q)
    c = mk.coalesce(b)
    mine = {"aligned_width_target_N": len(b["ts"]),
            "aligned_width_target_Total": (b["te"] - b["ts"] + 1).sum(),
            "aligned_width_query_Total": (b["qe"] - b["qs"] + 1).sum(),
            "chain_width_target_N": len(c["ts"]),
            "chain_width_target_Total": (c["te"] - c["ts"] + 1).sum(),
            "chain_width_query_Total": (c["qe"] - c["qs"] + 1).sum(),
            "guessed_target_length_N": len(np.unique(b["tseq"])),
            "guessed_query_length_N": len(np.unique(b["qseq"]))}
    for k in keys:
        ref = ex.iloc[i][k]
        d = abs(mine[k] - ref) / max(abs(ref), 1)
        worst[k] = max(worst[k], d)
        if d > 1e-6:
            bad.append((br, t[:25], q[:25], k, mine[k], ref))
print(worst)
print(len(bad), "discrepancies")
for x in bad[:20]:
    print(x)
