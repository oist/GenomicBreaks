"""Table of existing statistics (one row per ordered pair) from the YAML files
of the statistics pipeline, plus P recomputed from the count matrices."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import yaml
import mk_common as mk


def read_yaml(fn):
    with open(fn) as fh:
        d = yaml.safe_load(fh)
    if isinstance(d, str):
        d = yaml.safe_load(d)
    out = {}
    for k, v in d.items():
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            out[k] = float(v)
        else:
            try:
                out[k] = float(v)
            except (TypeError, ValueError):
                out[k] = np.nan
    return out


def p_matrix(fn):
    rows = [l.split() for l in open(fn) if l.strip() and not l.startswith("#")]
    cols = rows[0]
    same = tot = 0.0
    for r in rows[1:]:
        if r[0] in "ACGT":
            for letter, v in zip(cols, r[1:]):
                if letter in "ACGT":
                    tot += float(v)
                    if letter == r[0]:
                        same += float(v)
    return 1 - same / tot if tot else np.nan


if __name__ == "__main__":
    recs = []
    for br in mk.branches():
        for t, q in mk.pairs_of(br):
            d = read_yaml(mk.path(br, "yaml", t, q))
            d.update(branch=br, target=t, query=q,
                     P_matrix=p_matrix(mk.path(br, "matrix", t, q)))
            recs.append(d)
    df = pd.DataFrame(recs)
    front = ["branch", "target", "query"]
    df = df[front + [c for c in df.columns if c not in front]]
    df.to_csv(os.path.join(mk.RES, "existing_yaml.tsv"), sep="\t", index=False)
    print(df.shape)
    print(df.isna().mean().sort_values(ascending=False).head(40).to_string())
    print((df["P"] - df["P_matrix"]).abs().max())
