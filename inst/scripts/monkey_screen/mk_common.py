"""
Common code for the monkey-function screen: file layout, loading of
blasttabplus files, the 250-longest filter of the statistics pipeline,
coalescing of colinear blocks (re-implementation of
GenomicBreaks::coalesce_contigs with default arguments) and swap.

A pair is held in a dict of numpy arrays ("blocks"), sorted by target
sequence (order of first appearance) and target start:
  tseq, qseq   int32 codes of the sequences
  ts, te       target start and end (1-based, inclusive, ts <= te)
  qs, qe       query start and end
  strand       +1 / -1
  tL, qL       float64 lengths of all target / query sequences (by code)
and, for raw blocks only, the alignment columns of the file:
  pident, alen, mism, gapopen, bits, score
For coalesced blocks, the extra columns are:
  nmem (blocks merged), tcov, qcov (aligned bases in the chain on each
  genome), mism, gapopen, alen, score (sums), ident (matching bases).
"""
import gzip
import os
import numpy as np
import pandas as pd

ROOT = os.environ.get("MONKEY_ROOT", os.getcwd())
DATA = os.path.join(ROOT, "data", "monkey")
CACHE = os.path.join(ROOT, "cache")
RES = os.path.join(ROOT, "results")

COLS = ["qname", "tname", "pident", "alen", "mism", "gapopen",
        "qstart", "qend", "tstart", "tend", "E", "bits", "qlen", "tlen", "score"]
RAW_EXTRA = ["pident", "alen", "mism", "gapopen", "bits", "score"]
COA_EXTRA = ["nmem", "tcov", "qcov", "mism", "gapopen", "alen", "score", "ident"]


def branches():
    return sorted(b for b in os.listdir(DATA)
                  if os.path.isdir(os.path.join(DATA, b, "blasttabplus")))


def pairs_of(branch):
    d = os.path.join(DATA, branch, "blasttabplus")
    out = []
    for f in sorted(os.listdir(d)):
        if f.endswith(".blasttabplus.gz"):
            t, q = f[:-len(".blasttabplus.gz")].split("___")
            out.append((t, q))
    return out


def path(branch, kind, t, q):
    ext = {"blasttabplus": ".blasttabplus.gz", "matrix": ".o2o.matrix.txt",
           "o2o": ".o2o.tsv", "yaml": ".yaml"}[kind]
    return os.path.join(DATA, branch, kind, t + "___" + q + ext)


def read_blasttabplus(fn):
    """Raw blocks of a blasttabplus file, as load_genomic_breaks() does, then
    the filter of skeleton.Rmd (250 sequences with the most aligned bases)."""
    df = pd.read_csv(fn, sep="\t", header=None, names=COLS, compression="gzip",
                     dtype={"qname": str, "tname": str}, na_values=["NA"])
    tcode, tnames = pd.factorize(df["tname"], sort=False)
    qcode, qnames = pd.factorize(df["qname"], sort=False)
    tstart = df["tstart"].to_numpy(np.int64)
    tend = df["tend"].to_numpy(np.int64)
    qstart = df["qstart"].to_numpy(np.int64)
    qend = df["qend"].to_numpy(np.int64)
    assert (qstart <= qend).all(), "query on the minus strand in " + fn
    flip = tend < tstart
    b = dict(tseq=tcode.astype(np.int32), qseq=qcode.astype(np.int32),
             ts=np.where(flip, tend, tstart), te=np.where(flip, tstart, tend),
             qs=qstart, qe=qend, strand=np.where(flip, -1, 1).astype(np.int8))
    for c in RAW_EXTRA:
        b[c] = df[c].to_numpy(np.float64)
    tL = np.zeros(len(tnames))
    tL[tcode] = df["tlen"].to_numpy(np.float64)
    qL = np.zeros(len(qnames))
    qL[qcode] = df["qlen"].to_numpy(np.float64)
    b["tL"], b["qL"] = tL, qL
    b["tnames"], b["qnames"] = np.asarray(tnames, dtype=str), np.asarray(qnames, dtype=str)
    # --- filter of skeleton.Rmd
    if len(tnames) > 250 or len(qnames) > 250:
        tw = b["te"] - b["ts"] + 1
        s = np.bincount(b["tseq"], weights=tw, minlength=len(tnames))
        keep_t = np.argsort(-s, kind="stable")[:250]
        b = subset(b, np.isin(b["tseq"], keep_t))
        qw = b["qe"] - b["qs"] + 1
        s = np.bincount(b["qseq"], weights=qw, minlength=len(qnames))
        s[np.bincount(b["qseq"], minlength=len(qnames)) == 0] = -1
        keep_q = np.argsort(-s, kind="stable")[:250]
        keep_q = keep_q[s[keep_q] >= 0]
        b = subset(b, np.isin(b["qseq"], keep_q))
    return sort_target(b)


PER_BLOCK = ["tseq", "qseq", "ts", "te", "qs", "qe", "strand"]


def block_keys(b):
    return [k for k in b if k not in ("tL", "qL", "tnames", "qnames")]


def subset(b, mask):
    out = {k: (b[k][mask] if k in block_keys(b) else b[k]) for k in b}
    return out


def sort_target(b):
    o = np.lexsort((b["te"], b["ts"], b["tseq"]))
    return subset(b, o)


def swap(b):
    """Exchange the roles of target and query (GenomicBreaks::swap, then sort)."""
    out = dict(b)
    out["tseq"], out["qseq"] = b["qseq"], b["tseq"]
    out["ts"], out["te"], out["qs"], out["qe"] = b["qs"], b["qe"], b["ts"], b["te"]
    out["tL"], out["qL"] = b["qL"], b["tL"]
    if "tnames" in b:
        out["tnames"], out["qnames"] = b["qnames"], b["tnames"]
    if "tcov" in b:
        out["tcov"], out["qcov"] = b["qcov"], b["tcov"]
    return sort_target(out)


def query_order(b):
    """Permutation that sorts the blocks along the query genome."""
    return np.lexsort((b["qe"], b["qs"], b["qseq"]))


def colinear_with_next(b):
    """TRUE where a block (in target order) is colinear with the next one:
    same target sequence, same strand, same query sequence, and the next
    block is the immediate neighbour on the query in the direction given by
    the strand (flagColinearAlignments with tol = Inf)."""
    n = len(b["ts"])
    col = np.zeros(n, dtype=bool)
    if n < 2:
        return col
    qo = query_order(b)
    qrank = np.empty(n, dtype=np.int64)
    qrank[qo] = np.arange(n)
    same_t = b["tseq"][1:] == b["tseq"][:-1]
    same_q = b["qseq"][1:] == b["qseq"][:-1]
    same_s = b["strand"][1:] == b["strand"][:-1]
    step = qrank[1:] - qrank[:-1]
    col[:-1] = same_t & same_q & same_s & (step == b["strand"][:-1])
    return col


def coalesce(b):
    """Merge runs of colinear blocks.  Target and query ranges span the run."""
    n = len(b["ts"])
    col = colinear_with_next(b)
    start = np.ones(n, dtype=bool)
    start[1:] = ~col[:-1]
    idx = np.flatnonzero(start)
    out = dict(tL=b["tL"], qL=b["qL"])
    if "tnames" in b:
        out["tnames"], out["qnames"] = b["tnames"], b["qnames"]
    out["tseq"] = b["tseq"][idx]
    out["qseq"] = b["qseq"][idx]
    out["strand"] = b["strand"][idx]
    out["ts"] = np.minimum.reduceat(b["ts"], idx) if n else b["ts"]
    out["te"] = np.maximum.reduceat(b["te"], idx) if n else b["te"]
    out["qs"] = np.minimum.reduceat(b["qs"], idx) if n else b["qs"]
    out["qe"] = np.maximum.reduceat(b["qe"], idx) if n else b["qe"]
    tw = (b["te"] - b["ts"] + 1).astype(np.float64)
    qw = (b["qe"] - b["qs"] + 1).astype(np.float64)
    add = lambda x: np.add.reduceat(x, idx) if n else x
    out["nmem"] = add(np.ones(n))
    out["tcov"] = add(tw)
    out["qcov"] = add(qw)
    out["mism"] = add(b["mism"])
    out["gapopen"] = add(b["gapopen"])
    out["alen"] = add(b["alen"])
    out["score"] = add(b["score"])
    # identical bases, recovered exactly: aligned columns minus mismatches and gap columns
    out["ident"] = add(tw + qw - b["alen"] - b["mism"])
    return sort_target(out)


def cache_file(branch, t, q):
    return os.path.join(CACHE, branch, t + "___" + q + ".npz")


def save_cache(branch, t, q, b):
    os.makedirs(os.path.join(CACHE, branch), exist_ok=True)
    np.savez(cache_file(branch, t, q), **b)


def load_cache(branch, t, q):
    z = np.load(cache_file(branch, t, q), allow_pickle=False)
    return {k: z[k] for k in z.files}
