"""Reference table of existing statistics: YAML values, package indices that
are not in the YAML files, and simple quantities derived from them."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
import numpy as np, pandas as pd
import mk_common as mk, mk_existing as me


def one(args):
    br, t, q = args
    b = mk.load_cache(br, t, q)
    sw = mk.swap(b)
    c = mk.coalesce(b)
    bp, dist = me.bpgraph_totals(c)
    ng = len(np.unique(c["tseq"].astype(np.int64) * (int(c["qseq"].max()) + 1) + c["qseq"]))
    return dict(branch=br, target=t, query=q,
                index_karyotype_target=me.karyotype_index(b), index_karyotype_query=me.karyotype_index(sw),
                index_tau_target=me.tau_index(b), index_tau_query=me.tau_index(sw),
                bpg_breakpoints=bp, bpg_dist=dist,
                bpg_breakpoints_frac=bp / (len(b["ts"]) + ng), bpg_dist_frac=dist / (len(b["ts"]) + ng),
                bpg_reuse=2 * dist / bp if bp > 0 else np.nan)


if __name__ == "__main__":
    jobs = [(br, t, q) for br in mk.branches() for (t, q) in mk.pairs_of(br)]
    t0 = time.time()
    with Pool(2) as pool:
        recs = pool.map(one, jobs, chunksize=8)
    x = pd.DataFrame(recs)
    y = pd.read_csv(os.path.join(mk.RES, "existing_yaml.tsv"), sep="\t")
    d = y.merge(x, on=["branch", "target", "query"], validate="one_to_one")
    assert len(d) == len(y)
    # counts and totals of features that do not occur are zero, not missing
    for c in d.columns:
        if c.endswith("_N") or c.endswith("_Total"):
            d[c] = d[c].fillna(0)
    # quantities derived in ScrambledTreeBuilder::formatStats and the Branch notebooks
    d["der_percent_aligned"] = 100 * d["aligned_width_target_Total"] / d["guessed_target_length_Total"]
    d["der_avg_strandRand"] = (d["index_strandRand_target"] + d["index_strandRand_query"]) / 2
    d["der_avg_synteny"] = (d["index_synteny_target"] + d["index_synteny_query"]) / 2
    d["der_avg_chr_num"] = (d["guessed_target_length_N"] + d["guessed_query_length_N"]) / 2
    # contrasts between the target and the query versions of the same statistic
    lr = lambda a, b: np.log(d[a] / d[b])
    d["der_lr_length"] = lr("TargetLength", "QueryLength")
    d["der_lr_nseq"] = lr("guessed_target_length_N", "guessed_query_length_N")
    for s in ("Mean", "Median", "Max", "L50", "Total"):
        d["der_lr_seqlen_" + s] = lr("guessed_target_length_" + s, "guessed_query_length_" + s)
        d["der_lr_aligned_" + s] = lr("aligned_width_target_" + s, "aligned_width_query_" + s)
        d["der_lr_chain_" + s] = lr("chain_width_target_" + s, "chain_width_query_" + s)
    for s in ("synteny", "correlation", "GOCvicinity4", "strandRand", "karyotype", "tau"):
        d["der_d_" + s] = d["index_%s_target" % s] - d["index_%s_query" % s]
    d.to_csv(os.path.join(mk.RES, "existing_all.tsv"), sep="\t", index=False)
    print(d.shape, round(time.time() - t0), "s")
