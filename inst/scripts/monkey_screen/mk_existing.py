"""Python versions of GenomicBreaks functions that return one number per
alignment but are not in the YAML files of the statistics pipeline.  They
are used only as references to de-correlate against."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def karyotype_index(b):
    n = len(b["ts"])
    if n == 0:
        return 0.0
    if n == 1:
        return 0.0
    qw = (b["qe"] - b["qs"] + 1).astype(float)
    rows = []                                   # (target, query, value)
    for t in np.unique(b["tseq"]):
        m = b["tseq"] == t
        chr_length = b["te"][m].max() - b["ts"][m].min()
        qs_, inv = np.unique(b["qseq"][m], return_inverse=True)
        allq = np.bincount(inv, weights=qw[m])
        main = qs_[allq >= chr_length * 0.1]
        for q in main:
            rows.append((t, q, 1.0 / len(main)))
    if not rows:
        return 1.0
    rows = np.array(rows)
    vals = rows[:, 2].copy()
    for q in np.unique(rows[:, 1]):
        m = rows[:, 1] == q
        s = vals[m].sum()
        if s > 1:
            vals[m] = vals[m] / s
    return 1 - round(float(vals.mean()), 4)


def tau_index(b):
    n = len(b["ts"])
    if n == 1:
        return 1.0
    x = np.bincount(b["qseq"])
    x = x[x > 0].astype(float)
    if len(x) == 1:
        return 1.0
    return float(np.sum(1 - x / x.max()) / (len(x) - 1))


def bpgraph_totals(b):
    """Sums, over pairs of sequences, of the number of breakpoints and of
    N - cycles of the breakpoint graph (breakpointGraphProperties())."""
    n = len(b["ts"])
    nq = int(b["qseq"].max()) + 1
    gkey = b["tseq"].astype(np.int64) * nq + b["qseq"]
    gu, gid = np.unique(gkey, return_inverse=True)
    cnt = np.bincount(gid)
    start = np.r_[0, np.cumsum(cnt)[:-1]]
    o_t = np.lexsort((b["ts"], gid))
    o_q = np.lexsort((b["qe"], b["qs"], gid))
    qr = np.empty(n, dtype=np.int64)
    qr[o_q] = np.arange(n) - start[gid[o_q]] + 1        # 1-based rank on the query
    g = gid[o_t]
    pi = qr[o_t]
    sg = b["strand"][o_t]
    # extended permutation, per group, with a global offset per group
    off = (2 * start + 2 * np.arange(len(cnt)))          # offset of the values of each group
    left = np.where(sg > 0, 2 * pi, 2 * pi + 1)
    right = np.where(sg > 0, 2 * pi + 1, 2 * pi)
    # sequence of values: 1, left_1, right_1, ..., left_n, right_n, 2n+2
    first = np.r_[True, g[1:] != g[:-1]]
    last = np.r_[g[1:] != g[:-1], True]
    prev_right = np.r_[0, right[:-1]]
    prev_right[first] = 1
    a = prev_right                     # black edge (a, left) before each block
    bl = left
    end_a = right[last]
    end_b = 2 * cnt[g[last]] + 2
    A = np.r_[a + off[g] - 1, end_a + off[g[last]] - 1]
    B = np.r_[bl + off[g] - 1, end_b + off[g[last]] - 1]
    bp = np.bincount(np.r_[g, g[last]], weights=(np.abs(np.r_[a, end_a] - np.r_[bl, end_b]) != 1),
                     minlength=len(cnt))
    # gray edges: values (2i+1, 2i+2) for i = 0..n
    tot = int(2 * n + 2 * len(cnt))
    gs = np.repeat(off, cnt + 1) + 2 * (np.arange(n + len(cnt)) - np.repeat(start + np.arange(len(cnt)), cnt + 1))
    G1, G2 = gs, gs + 1
    M = coo_matrix((np.ones(len(A) + len(G1), dtype=np.int8), (np.r_[A, G1], np.r_[B, G2])), shape=(tot, tot))
    ncomp, lab = connected_components(M, directed=False)
    firstv = np.unique(lab, return_index=True)[1]
    grp_of_vertex = np.repeat(np.arange(len(cnt)), 2 * cnt + 2)
    cyc = np.bincount(grp_of_vertex[firstv], minlength=len(cnt))
    N = cnt + 1
    return float(bp.sum()), float((N - cyc).sum())
