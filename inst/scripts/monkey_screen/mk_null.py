"""Transformations of an alignment used by the screen:
 - shuffle(): random arrangement (order along the query and strands), the
   blocks keeping their target position, widths and alignment columns;
 - rerepresent(): the same alignment written differently (sequences reversed
   at random, sequence codes permuted).  A sound function gives the same
   value before and after."""
import numpy as np
import mk_common as mk


def _relayout(seq, start, end, L, width, rng):
    """Random order of the blocks along one genome.  The genome keeps its
    slots (number of blocks per sequence and the gap before each slot); the
    blocks, with their widths, are dealt at random into the slots.
    Returns new sequence codes, starts, ends (per block) and lengths."""
    n = len(seq)
    o = np.lexsort((end, start, seq))
    seq_s, s_s, e_s = seq[o], start[o], end[o]
    first = np.r_[True, seq_s[1:] != seq_s[:-1]]
    prev_end = np.r_[0, e_s[:-1]]
    prev_end[first] = 0
    gap_before = s_s - prev_end - 1
    perm = rng.permutation(n)                 # slot k receives block perm[k]
    w = width[perm]
    cs = np.cumsum(gap_before + w)
    idx = np.flatnonzero(first)
    offset = np.repeat(np.r_[0, cs[idx[1:] - 1]], np.diff(np.r_[idx, n]))
    new_e = cs - offset
    new_s = new_e - w + 1
    nseq = np.empty_like(seq); ns = np.empty_like(start); ne = np.empty_like(end)
    nseq[perm], ns[perm], ne[perm] = seq_s, new_s, new_e
    last = np.r_[idx[1:] - 1, n - 1]
    newL = np.asarray(L, dtype=float).copy()
    newL[seq_s[last]] = new_e[last] + (L[seq_s[last]] - e_s[last])
    return nseq, ns, ne, newL


def shuffle_both(b, rng):
    """Null arrangement: the order of the blocks is random along the target
    and, independently, along the query, and the strands are dealt at random
    among the blocks.  Each block keeps its widths and alignment columns."""
    out = dict(b)
    n = len(b["ts"])
    if n < 2:
        return out
    tw = b["te"] - b["ts"] + 1
    qw = b["qe"] - b["qs"] + 1
    out["tseq"], out["ts"], out["te"], out["tL"] = _relayout(b["tseq"], b["ts"], b["te"], b["tL"], tw, rng)
    out["qseq"], out["qs"], out["qe"], out["qL"] = _relayout(b["qseq"], b["qs"], b["qe"], b["qL"], qw, rng)
    out["strand"] = rng.permutation(b["strand"])
    return mk.sort_target(out)


def shuffle(b, rng):
    """As shuffle_both(), but the target side is left untouched."""
    n = len(b["ts"])
    out = dict(b)
    if n < 2:
        return out
    qo = mk.query_order(b)
    qseq_s, qs_s, qe_s = b["qseq"][qo], b["qs"][qo], b["qe"][qo]
    first = np.r_[True, qseq_s[1:] != qseq_s[:-1]]
    prev_end = np.r_[0, qe_s[:-1]]
    prev_end[first] = 0
    gap_before = qs_s - prev_end - 1
    perm = rng.permutation(n)                 # slot k receives block perm[k]
    w = (b["qe"] - b["qs"] + 1)[perm]
    cs = np.cumsum(gap_before + w)
    idx = np.flatnonzero(first)
    sizes = np.diff(np.r_[idx, n])
    offset = np.repeat(np.r_[0, cs[idx[1:] - 1]], sizes)
    new_qe = cs - offset
    new_qs = new_qe - w + 1
    qseq = np.empty_like(b["qseq"]); qs = np.empty_like(b["qs"]); qe = np.empty_like(b["qe"])
    qseq[perm], qs[perm], qe[perm] = qseq_s, new_qs, new_qe
    out["qseq"], out["qs"], out["qe"] = qseq, qs, qe
    # sequence lengths: keep the trailing unaligned part of each query sequence
    last = np.r_[idx[1:] - 1, n - 1]
    qL = b["qL"].astype(float).copy()
    qL[qseq_s[last]] = new_qe[last] + (b["qL"][qseq_s[last]] - qe_s[last])
    out["qL"] = qL
    out["strand"] = rng.permutation(b["strand"])
    return out


def rerepresent(b, rng):
    out = dict(b)
    nT, nQ = len(b["tL"]), len(b["qL"])
    ft = rng.random(nT) < 0.5
    fq = rng.random(nQ) < 0.5
    s = b["strand"].copy()
    ts, te, qs, qe = (b[k].astype(np.int64).copy() for k in ("ts", "te", "qs", "qe"))
    m = ft[b["tseq"]]
    L = b["tL"][b["tseq"]].astype(np.int64)
    ts[m], te[m] = (L - b["te"] + 1)[m], (L - b["ts"] + 1)[m]
    s[m] = -s[m]
    m = fq[b["qseq"]]
    L = b["qL"][b["qseq"]].astype(np.int64)
    qs[m], qe[m] = (L - b["qe"] + 1)[m], (L - b["qs"] + 1)[m]
    s[m] = -s[m]
    pt, pq = rng.permutation(nT), rng.permutation(nQ)
    out["tseq"] = pt[b["tseq"]].astype(b["tseq"].dtype)
    out["qseq"] = pq[b["qseq"]].astype(b["qseq"].dtype)
    tL = np.empty(nT); tL[pt] = b["tL"]
    qL = np.empty(nQ); qL[pq] = b["qL"]
    out.update(ts=ts, te=te, qs=qs, qe=qe, strand=s, tL=tL, qL=qL)
    for k in ("tnames", "qnames"):
        out.pop(k, None)
    return mk.sort_target(out)
