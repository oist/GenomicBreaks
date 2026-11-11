"""Screening statistics, vectorised over candidates.

Layout: one row per ordered pair (Branch, target, query), in a fixed order
given by `load_design()`; candidate values in matrices rows x candidates."""
import os
import numpy as np
import pandas as pd
from scipy.stats import rankdata
import mk_common as mk


# ------------------------------------------------------------------ design

def upgma_nodes(ids, d):
    """Average-linkage clustering of `ids` with distances d[(a, b)].
    Returns {frozenset((a, b)): node index}; the node of a pair is the merge
    that first joins the two genomes (their last common ancestor)."""
    clusters = {i: [g] for i, g in enumerate(ids)}
    dist = {}
    keys = list(clusters)
    for x, i in enumerate(keys):
        for j in keys[x + 1:]:
            dist[(i, j)] = d[frozenset((ids[i], ids[j]))]
    node, nxt, k = {}, len(ids), 0
    heights = []
    while len(clusters) > 1:
        (i, j), h = min(dist.items(), key=lambda kv: (kv[1], kv[0]))
        for a in clusters[i]:
            for b in clusters[j]:
                node[frozenset((a, b))] = k
        heights.append(h)
        ni, nj = len(clusters[i]), len(clusters[j])
        new = {}
        for c in clusters:
            if c in (i, j):
                continue
            new[c] = (dist[(min(i, c), max(i, c))] * ni + dist[(min(j, c), max(j, c))] * nj) / (ni + nj)
        for c in list(dist):
            if i in c or j in c:
                del dist[c]
        members = clusters.pop(i) + clusters.pop(j)
        for c, v in new.items():
            dist[(c, nxt)] = v
        clusters[nxt] = members
        nxt += 1
        k += 1
    return node, heights


def load_design():
    """Table of ordered pairs with: mate (row of the reciprocal pair), upair
    (index of the unordered pair), node (ancestor node, unique across
    Branches), P (mean of both orders)."""
    ex = pd.read_csv(os.path.join(mk.RES, "existing_all.tsv"), sep="\t")
    key = {(r.branch, r.target, r.query): i for i, r in enumerate(ex.itertuples())}
    ex["mate"] = [key[(r.branch, r.query, r.target)] for r in ex.itertuples()]
    ex["Pm"] = (ex["P"].to_numpy() + ex["P"].to_numpy()[ex["mate"].to_numpy()]) / 2
    up, node, gi = {}, np.zeros(len(ex), dtype=int), {}
    upair = np.zeros(len(ex), dtype=int)
    off = 0
    for br, d in ex.groupby("branch", sort=False):
        ids = sorted(set(d.target))
        dd = {frozenset((r.target, r.query)): r.Pm for r in d.itertuples()}
        nd, _ = upgma_nodes(ids, dd)
        for r in d.itertuples():
            fs = frozenset((r.target, r.query))
            node[r.Index] = off + nd[fs]
            upair[r.Index] = up.setdefault((br, fs), len(up))
        off += len(ids) - 1
    ex["node"] = node
    ex["upair"] = upair
    return ex


def load_stage1(design, sub="stage1"):
    d = os.path.join(mk.RES, sub)
    names = [l.strip() for l in open(os.path.join(d, "names.txt"))]
    fwd = np.full((len(design), len(names)), np.nan)
    swp = np.full((len(design), len(names)), np.nan)
    key = {(r.branch, r.target, r.query): i for i, r in enumerate(design.itertuples())}
    for br in design.branch.unique():
        z = np.load(os.path.join(d, br + ".npz"))
        idx = [key[(br, t, q)] for t, q in zip(z["target"], z["query"])]
        fwd[idx] = z["fwd"]
        swp[idx] = z["swp"]
    return names, fwd, swp


# ----------------------------------------------------------------- helpers

def rank_cols(X):
    """Column-wise average ranks; NaN replaced by the mean rank of the column."""
    R = np.empty_like(X, dtype=float)
    for j in range(X.shape[1]):
        x = X[:, j]
        fin = np.isfinite(x)
        r = np.full(x.size, np.nan)
        if fin.any():
            r[fin] = rankdata(x[fin])
            r[~fin] = r[fin].mean()
        R[:, j] = r
    return R


def zs(R):
    """Centre and scale columns to unit sum of squares; constant columns -> 0."""
    R = R - np.nanmean(R, axis=0, keepdims=True)
    ss = np.sqrt(np.nansum(R * R, axis=0, keepdims=True))
    with np.errstate(invalid="ignore", divide="ignore"):
        Z = R / ss
    Z[~np.isfinite(Z)] = 0.0
    return Z


def colcorr(A, B):
    """Pearson correlation between matching columns of A and B (NaN-aware)."""
    ok = np.isfinite(A) & np.isfinite(B)
    A = np.where(ok, A, np.nan)
    B = np.where(ok, B, np.nan)
    A = A - np.nanmean(A, axis=0)
    B = B - np.nanmean(B, axis=0)
    num = np.nansum(A * B, axis=0)
    den = np.sqrt(np.nansum(A * A, axis=0) * np.nansum(B * B, axis=0))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = num / den
    r[ok.sum(axis=0) < 8] = np.nan
    return r


def spearman_cols(A, B):
    """Spearman correlation between matching columns, pairwise complete."""
    ok = np.isfinite(A) & np.isfinite(B)
    RA = np.full(A.shape, np.nan)
    RB = np.full(B.shape, np.nan)
    for j in range(A.shape[1]):
        m = ok[:, j]
        if m.sum() >= 8:
            RA[m, j] = rankdata(A[m, j])
            RB[m, j] = rankdata(B[m, j])
    return colcorr(RA, RB)


def partial(rxy, rxz, ryz, limit=0.98):
    """Partial correlation of x and y given z.  When x or y is itself almost
    a function of z (|r| above `limit`), nothing is left to correlate and the
    ratio is numerically unstable: the result is set to 0."""
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (rxy - rxz * ryz) / np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    bad = (np.abs(rxz) > limit) | (np.abs(ryz) > limit) | ~np.isfinite(r)
    return np.where(bad, 0.0, r)


# ---------------------------------------------------------- screen metrics

def distribution_metrics(V, branch):
    fin = np.isfinite(V)
    out = {"fin": fin.mean(axis=0)}
    out["fin_min"] = np.min([fin[branch == b].mean(axis=0) for b in np.unique(branch)], axis=0)
    nd = np.zeros(V.shape[1])
    mode = np.ones(V.shape[1])
    for j in range(V.shape[1]):
        x = V[fin[:, j], j]
        if x.size:
            with np.errstate(divide="ignore", invalid="ignore"):
                mag = np.where(x != 0, np.floor(np.log10(np.abs(x))), 0)
                y = np.round(x / 10.0 ** mag, 8) * 10.0 ** mag     # 9 significant digits
            u, c = np.unique(y, return_counts=True)
            nd[j] = u.size / x.size
            mode[j] = c.max() / x.size
    out["distinct"] = nd
    out["mode"] = mode
    return out


def reciprocal_metrics(fwd, swp, design):
    """Agreement between f(T___Q) and f(swap(Q___T)), which describe the same
    genome pair with the same roles but come from two alignments; and
    symmetry, between f(T___Q) and f(Q___T)."""
    mate = design["mate"].to_numpy()
    branch = design["branch"].to_numpy()
    rec, sym = [], []
    for b in pd.unique(branch):
        m = branch == b
        rec.append(spearman_cols(fwd[m], swp[mate][m]))
        sym.append(spearman_cols(fwd[m], fwd[mate][m]))
    rec, sym = np.array(rec), np.array(sym)
    with np.errstate(all="ignore"):
        out = {"recip_med": np.nanmedian(rec, axis=0), "recip_min": np.nanmin(rec, axis=0),
               "recip_q25": np.nanquantile(rec, 0.25, axis=0),
               "recip_pool": spearman_cols(fwd, swp[mate]),
               "sym_med": np.nanmedian(sym, axis=0), "sym_pool": spearman_cols(fwd, fwd[mate])}
    return out


def unordered(V, design):
    """Mean of both orders, one row per unordered pair; with Branch, node, P."""
    up = design["upair"].to_numpy()
    first = np.unique(up, return_index=True)[1]
    mate = design["mate"].to_numpy()
    U = (V[first] + V[mate[first]]) / 2
    d = design.iloc[first]
    return U, d["branch"].to_numpy(), d["node"].to_numpy(), d["Pm"].to_numpy(), d


def node_metrics(V, design, nperm=300, seed=1):
    """Similarity of the pairs that share an ancestor node, per Branch:
    share of the variance of the ranks between nodes (eta2), compared with
    permutations of the pairs within the Branch (z).  Same after removing
    the linear effect of the rank of P ("R")."""
    U, ubr, unode, uP, ud = unordered(V, design)
    rng = np.random.default_rng(seed)
    res = {k: [] for k in ("eta", "z", "etaR", "zR", "gen", "zgen")}
    targets, queries = ud["target"].to_numpy(), ud["query"].to_numpy()
    for b in pd.unique(ubr):
        m = ubr == b
        R = rank_cols(U[m])
        n = R.shape[0]
        nodes, inv = np.unique(unode[m], return_inverse=True)
        Z = np.zeros((n, nodes.size))
        Z[np.arange(n), inv] = 1
        cnt = Z.sum(0)
        # genome incidence (additive effect of each genome)
        gs = sorted(set(targets[m]) | set(queries[m]))
        gi = {g: k for k, g in enumerate(gs)}
        G = np.zeros((n, len(gs)))
        G[np.arange(n), [gi[t] for t in targets[m]]] = 1
        G[np.arange(n), [gi[q] for q in queries[m]]] = 1
        G = G - G.mean(0)
        Qg = np.linalg.qr(G)[0][:, :np.linalg.matrix_rank(G)]
        rp = rankdata(uP[m])
        rp = (rp - rp.mean()) / np.sqrt(np.sum((rp - rp.mean()) ** 2))

        def stats(M):
            M = M - M.mean(0)
            sst = np.sum(M * M, axis=0)
            with np.errstate(invalid="ignore", divide="ignore"):
                eta = np.sum((Z.T @ M) ** 2 / cnt[:, None], axis=0) / sst
                gen = np.sum((Qg.T @ M) ** 2, axis=0) / sst
            return eta, gen

        Rc = R - R.mean(0)
        RR = Rc - np.outer(rp, rp @ Rc)                 # residual on the rank of P
        e0, g0 = stats(R)
        e0R, _ = stats(RR)
        pe, peR, pg = [], [], []
        for _ in range(nperm):
            p = rng.permutation(n)
            a, c = stats(R[p])
            pe.append(a)
            pg.append(c)
            peR.append(stats(RR[p])[0])
        pe, peR, pg = np.array(pe), np.array(peR), np.array(pg)
        with np.errstate(invalid="ignore", divide="ignore"):
            res["eta"].append(e0)
            res["z"].append((e0 - pe.mean(0)) / pe.std(0))
            res["etaR"].append(e0R)
            res["zR"].append((e0R - peR.mean(0)) / peR.std(0))
            res["gen"].append(g0)
            res["zgen"].append((g0 - pg.mean(0)) / pg.std(0))
    out = {}
    nb = len(res["eta"])
    with np.errstate(all="ignore"):
        out["node_eta"] = np.nanmedian(res["eta"], axis=0)
        out["node_z"] = np.nansum(res["z"], axis=0) / np.sqrt(nb)
        out["nodeR_eta"] = np.nanmedian(res["etaR"], axis=0)
        out["nodeR_z"] = np.nansum(res["zR"], axis=0) / np.sqrt(nb)
        out["genome_r2"] = np.nanmedian(res["gen"], axis=0)
        out["genome_z"] = np.nansum(res["zgen"], axis=0) / np.sqrt(nb)
    return out


def corr_blocks(V, E, design, pcol="P"):
    """Rank correlations between candidates (columns of V) and reference
    statistics (columns of E): pooled, pooled partial on P, median of the
    within-Branch values, and the same partial on P.
    Returns four arrays (candidates x references)."""
    branch = design["branch"].to_numpy()
    P = design[pcol].to_numpy()

    def block(rows):
        RV = zs(rank_cols(V[rows]))
        rp = zs(rank_cols(P[rows][:, None]))[:, 0]
        cvp = RV.T @ rp
        C = np.full((V.shape[1], E.shape[1]), np.nan)
        Cp = np.full((V.shape[1], E.shape[1]), np.nan)
        Efin = np.isfinite(E[rows])
        full = Efin.all(axis=0)
        if full.any():
            RE = zs(rank_cols(E[rows][:, full]))
            c = RV.T @ RE
            C[:, full] = c
            Cp[:, full] = partial(c, cvp[:, None], (RE.T @ rp)[None, :])
        for j in np.flatnonzero(~full):
            m = Efin[:, j]
            if m.sum() < max(12, 0.2 * m.size):
                continue
            rv = zs(rank_cols(V[rows][m]))
            re = zs(rank_cols(E[rows][m, j:j + 1]))[:, 0]
            rpp = zs(rank_cols(P[rows][m][:, None]))[:, 0]
            c = rv.T @ re
            C[:, j] = c
            Cp[:, j] = partial(c, rv.T @ rpp, re @ rpp)
        return C, Cp

    pooled, pooled_p = block(np.ones(len(branch), bool))
    w, wp = [], []
    for b in pd.unique(branch):
        c, cp = block(branch == b)
        w.append(c)
        wp.append(cp)
    with np.errstate(all="ignore"):
        within = np.nanmedian(np.array(w), axis=0)
        within_p = np.nanmedian(np.array(wp), axis=0)
    return pooled, pooled_p, within, within_p
