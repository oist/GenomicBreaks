"""
Monkey-typed summary functions of a pairwise genome alignment.

A candidate function is the combination of
  - an input level: raw one-to-one blocks ("raw") or coalesced blocks ("coa");
  - one or two per-block or per-junction vectors;
  - a reducer that turns them into one number.
Nothing here encodes a biological hypothesis: the vectors are the columns of
the alignment table and elementary combinations of them, and the reducers
are generic descriptive statistics.  No candidate takes a parameter.

Vector classes
  B   one value per block (target order)
  J   one value per junction: two consecutive blocks on the same target
      sequence
  JQ  junctions whose two blocks are also on the same query sequence

All vectors are built so that the result does not depend on arbitrary
choices of representation: the order of the sequences in the files and the
orientation (forward or reverse-complement) of each sequence.

Candidate names:  <level>.<reducer>.<vector>[~<vector>]
"""
import numpy as np
from bisect import bisect_left
from scipy.stats import rankdata
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

NAN = float("nan")
TOL = 1e-18      # a variance below TOL * (mean^2 + 1) counts as zero


def _zero(ss, n, m):
    """TRUE when the sum of squares `ss` of `n` values of mean `m` is only
    floating-point noise."""
    return ss <= TOL * n * (m * m + 1.0)


# ------------------------------------------------------------------ helpers

def _moments(y):
    n = y.size
    if n < 3:
        return NAN, NAN, NAN, NAN
    m = y.mean()
    d = y - m
    v = np.dot(d, d) / n
    if _zero(v * n, n, m):
        return m, 0.0, NAN, NAN
    sd = np.sqrt(v)
    d2 = d * d
    return m, sd, np.dot(d2, d) / n / sd ** 3, np.dot(d2, d2) / n / v ** 2 - 3.0


def _q(s, p):
    """Quantile of a sorted array, linear interpolation (R type 7)."""
    h = (s.size - 1) * p
    lo = int(np.floor(h))
    hi = min(lo + 1, s.size - 1)
    return s[lo] + (h - lo) * (s[hi] - s[lo])


def _pearson(a, b):
    n = a.size
    if n < 3:
        return NAN
    ma, mb = a.mean(), b.mean()
    a = a - ma
    b = b - mb
    sa, sb = np.dot(a, a), np.dot(b, b)
    if _zero(sa, n, ma) or _zero(sb, n, mb):
        return NAN
    return float(np.dot(a, b) / np.sqrt(sa * sb))


def _eta2(y, g):
    """Share of the variance of y that lies between the groups g."""
    if y.size < 3:
        return NAN
    m = y.mean()
    sst = np.dot(y - m, y - m)
    if _zero(sst, y.size, m):
        return NAN
    _, inv = np.unique(g, return_inverse=True)
    cnt = np.bincount(inv)
    if cnt.size < 2:
        return 0.0
    mg = np.bincount(inv, weights=y) / cnt
    return float(np.dot(cnt, (mg - m) ** 2) / sst)


def _acf(y, grp, lag):
    """Correlation between values `lag` positions apart on the same sequence.
    The two members of a pair are exchangeable (one mean and one variance for
    both), so that the result does not depend on the reading direction."""
    if y.size <= lag + 2:
        return NAN
    ok = grp[lag:] == grp[:-lag]
    a, b = y[:-lag][ok], y[lag:][ok]
    fin = np.isfinite(a) & np.isfinite(b)
    if not fin.all():
        a, b = a[fin], b[fin]
    if a.size < 3:
        return NAN
    m = (a.sum() + b.sum()) / (2 * a.size)
    a = a - m
    b = b - m
    den = (np.dot(a, a) + np.dot(b, b)) / 2
    if _zero(den, a.size, m):
        return NAN
    return float(np.dot(a, b) / den)


def _turn(y, grp):
    """Proportion of interior positions that are strict local extrema."""
    if y.size < 5:
        return NAN
    ok = (grp[2:] == grp[:-2]) & np.isfinite(y[2:]) & np.isfinite(y[1:-1]) & np.isfinite(y[:-2])
    if ok.sum() < 3:
        return NAN
    d1 = y[1:-1] - y[:-2]
    d2 = y[2:] - y[1:-1]
    return float((((d1 * d2) < 0) & ok).sum() / ok.sum())


# ------------------------------------------------------------------ context

class Ctx:
    """Vectors of one alignment at one level, computed on demand."""

    CORE = ["tw", "ident", "sdens", "asym", "tedge", "qedge"]

    def __init__(self, b, level):
        self.level = level
        self.b = b
        self.n = n = len(b["ts"])
        f = np.float64
        self.ts, self.te = b["ts"].astype(f), b["te"].astype(f)
        self.qs, self.qe = b["qs"].astype(f), b["qe"].astype(f)
        self.tseq, self.qseq = b["tseq"].astype(np.int64), b["qseq"].astype(np.int64)
        self.strand = b["strand"].astype(f)
        self.tw = self.te - self.ts + 1
        self.qw = self.qe - self.qs + 1
        self._v = {}
        # query order and ranks
        qo = np.lexsort((b["qe"], b["qs"], b["qseq"]))
        self.qorder = qo
        self.qrank = np.empty(n, dtype=np.int64)
        self.qrank[qo] = np.arange(n)
        # pairs of sequences ("groups")
        nq = int(self.qseq.max()) + 1 if n else 1
        gkey = self.tseq * nq + self.qseq
        self.gu, self.gid = np.unique(gkey, return_inverse=True)
        self.gW = np.bincount(self.gid, weights=self.tw)
        # Majority strand of each pair of sequences, by aligned target bases.
        # Exact ties happen (for instance two blocks of the same width on
        # opposite strands); they are broken with other block properties that
        # do not depend on the orientation of the sequences, so that the
        # answer is the same however the sequences are written.
        self.gmaj = np.ones(self.gu.size)
        todo = np.ones(self.gu.size, dtype=bool)
        for wname in ("tw", "qw", "tedge", "qedge"):
            if not todo.any():
                break
            w = self.tw if wname == "tw" else self.qw if wname == "qw" else None
            if w is None:
                L = b["tL"][self.tseq] if wname == "tedge" else b["qL"][self.qseq]
                mid = (self.ts + self.te) / 2 if wname == "tedge" else (self.qs + self.qe) / 2
                w = np.round(1 - np.abs(2 * (mid - 0.5) / L - 1), 12)
            gS = np.bincount(self.gid, weights=w * self.strand, minlength=self.gu.size)
            self.gmaj[todo & (gS < 0)] = -1.0
            todo &= gS == 0
        # junctions
        if n > 1:
            same_t = self.tseq[1:] == self.tseq[:-1]
            self.iJ = np.flatnonzero(same_t)
        else:
            self.iJ = np.zeros(0, dtype=np.int64)
        i = self.iJ
        self.jseq = self.tseq[i]
        self.jsameq = self.qseq[i] == self.qseq[i + 1] if i.size else np.zeros(0, bool)
        self.iJQ = i[self.jsameq]
        self.classes = {}
        self.kinds = {}

    # -- vector catalogue ------------------------------------------------
    def names(self, cls):
        if cls == "B":
            v = ["tw", "qw", "asym", "aasym", "tedge", "qedge", "tLb", "qLb", "lrat",
                 "smaj", "pdisp", "rdisp", "ident", "mismf", "gapf", "sdens", "gcol",
                 "bgap", "bflip", "bsameq", "tshare", "qshare"]
            v += ["alen", "glen"] if self.level == "raw" else ["nmem", "tdens", "qdens"]
            return v
        if cls == "J":
            core = self.CORE + (["nmem", "tdens"] if self.level == "coa" else [])
            return ["tgap", "flip", "sameq"] + ["jm_" + c for c in core] + ["jd_" + c for c in core]
        if cls == "JQ":
            return ["qstep", "rstep", "dirc", "gdiff", "colin"]

    def cls_of(self, name):
        for c in ("B", "J", "JQ"):
            if name in self.names(c):
                return c
        raise KeyError(name)

    def kind(self, name):
        if name in ("tpos", "qpos", "trk", "qrk"):
            return "aux"
        if name in ("tw", "qw", "tLb", "qLb", "alen", "glen", "nmem", "bgap", "tgap", "qstep", "rstep"):
            return "size"
        if name in ("smaj", "flip", "sameq", "colin"):
            return "bin"
        if name in ("bflip", "bsameq", "dirc"):
            return "few"
        return "frac"

    def canon(self, name):
        """Vector on its canonical scale: log for sizes, identity otherwise."""
        key = "c:" + name
        if key not in self._v:
            x = self.vec(name)
            self._v[key] = np.log(x) if self.kind(name) == "size" else x
        return self._v[key]

    def rank(self, name):
        """Average ranks among the finite values of the vector (NaN elsewhere)."""
        key = "r:" + name
        if key not in self._v:
            x = self.vec(name)
            fin = np.isfinite(x)
            if fin.all():
                r = rankdata(x)
            else:
                r = np.full(x.size, np.nan)
                if fin.any():
                    r[fin] = rankdata(x[fin])
            self._v[key] = r
        return self._v[key]

    def vec(self, name):
        if name not in self._v:
            x = np.asarray(self._make(name), dtype=np.float64)
            if self.kind(name) in ("frac", "few"):
                # ratios and positions: drop floating-point noise, which would
                # otherwise break ties differently in equivalent representations
                x = np.round(x, 12)
            self._v[name] = x
        return self._v[name]

    def _make(self, name):
        b, n = self.b, self.n
        tw, qw = self.tw, self.qw
        i = self.iJ
        if name == "tw":
            return tw
        if name == "qw":
            return qw
        if name == "asym":
            return (tw - qw) / (tw + qw)
        if name == "aasym":
            return np.abs(tw - qw) / (tw + qw)
        if name in ("tpos", "qpos"):
            # centre of the block, from 0 (first base) to 1 (last base)
            if name == "tpos":
                return ((self.ts + self.te) / 2 - 0.5) / b["tL"][self.tseq]
            return ((self.qs + self.qe) / 2 - 0.5) / b["qL"][self.qseq]
        if name == "tedge":
            return 1 - np.abs(2 * self.vec("tpos") - 1)
        if name == "qedge":
            return 1 - np.abs(2 * self.vec("qpos") - 1)
        if name == "tLb":
            return b["tL"][self.tseq]
        if name == "qLb":
            return b["qL"][self.qseq]
        if name == "lrat":
            return np.log(b["tL"][self.tseq]) - np.log(b["qL"][self.qseq])
        if name == "smaj":
            return self.strand * self.gmaj[self.gid]
        if name == "pdisp":
            qp = self.vec("qpos")
            qp = np.where(self.gmaj[self.gid] > 0, qp, 1 - qp)
            return np.abs(self.vec("tpos") - qp)
        if name in ("trk", "qrk"):
            if name == "trk":
                seq, order = self.tseq, np.arange(n)
            else:
                seq, order = self.qseq, self.qorder
            s = seq[order]
            first = np.r_[True, s[1:] != s[:-1]]
            start = np.maximum.accumulate(np.where(first, np.arange(n), 0))
            cnt = np.bincount(s)[s]
            r = np.empty(n)
            r[order] = (np.arange(n) - start + 0.5) / cnt
            return r
        if name == "rdisp":
            qr = self.vec("qrk")
            qr = np.where(self.gmaj[self.gid] > 0, qr, 1 - qr)
            return np.abs(self.vec("trk") - qr)
        if name == "tshare":
            return self.gW[self.gid] / np.bincount(self.tseq, weights=tw)[self.tseq]
        if name == "qshare":
            return self.gW[self.gid] / np.bincount(self.qseq, weights=tw)[self.qseq]
        # ---- alignment columns
        if self.level == "raw":
            alen, mism, gapo, score = b["alen"], b["mism"], b["gapopen"], b["score"]
            match = tw + qw - alen - mism
            gapc = 2 * alen - tw - qw
        else:
            alen, mism, gapo, score = b["alen"], b["mism"], b["gapopen"], b["score"]
            match = b["ident"]
            gapc = 2 * alen - b["tcov"] - b["qcov"]
        if name == "alen":
            return alen
        if name == "ident":
            return match / alen
        if name == "mismf":
            return mism / alen
        if name == "gapf":
            return gapo / alen
        if name == "sdens":
            return score / alen
        if name == "gcol":
            return gapc / alen
        if name == "glen":
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.where(gapo > 0, gapc / gapo, np.nan)
        if name == "nmem":
            return b["nmem"]
        if name == "tdens":
            return b["tcov"] / tw
        if name == "qdens":
            return b["qcov"] / qw
        # ---- junctions (compressed on the J set)
        if name == "tgap":
            return np.maximum(self.ts[i + 1] - self.te[i] - 1, 0) + 1
        if name == "flip":
            # strands are compared relative to the majority strand of each pair
            # of sequences, because the orientation of a sequence is arbitrary
            sm = self.vec("smaj")
            return (sm[i] != sm[i + 1]).astype(float)
        if name == "sameq":
            return self.jsameq.astype(float)
        if name.startswith("jm_"):
            y = self.canon(name[3:])
            return (y[i] + y[i + 1]) / 2
        if name.startswith("jd_"):
            y = self.canon(name[3:])
            return np.abs(y[i] - y[i + 1])
        # ---- blocks described by their flanking junctions
        if name in ("bgap", "bflip", "bsameq"):
            src = {"bgap": np.log(self.vec("tgap")), "bflip": self.vec("flip"),
                   "bsameq": self.vec("sameq")}[name]
            s = np.zeros(n)
            c = np.zeros(n)
            s[i] += src
            c[i] += 1
            s[i + 1] += src
            c[i + 1] += 1
            with np.errstate(invalid="ignore", divide="ignore"):
                m = np.where(c > 0, s / c, np.nan)
            return np.exp(m) if name == "bgap" else m
        # ---- junctions on the same query sequence (compressed on the JQ set)
        k = self.iJQ
        if name == "qstep":
            return np.maximum(np.maximum(self.qs[k], self.qs[k + 1])
                              - np.minimum(self.qe[k], self.qe[k + 1]) - 1, 0) + 1
        if name == "rstep":
            return np.abs(self.qrank[k + 1] - self.qrank[k]).astype(float)
        if name == "dirc":
            return np.sign(self.qrank[k + 1] - self.qrank[k]) * (self.strand[k] + self.strand[k + 1]) / 2
        if name == "gdiff":
            tg = np.maximum(self.ts[k + 1] - self.te[k] - 1, 0)
            qg = self.vec("qstep") - 1
            return (tg - qg) / (tg + qg + 2)
        if name == "colin":
            return ((self.vec("rstep") == 1) & (self.vec("dirc") == 1)).astype(float)
        raise KeyError(name)

    def grp(self, cls):
        return self.tseq if cls == "B" else self.jseq

    def on_jq(self, arr):
        """Restrict a J-class array to the JQ junctions."""
        return arr[self.jsameq]


# --------------------------------------------------------- single reducers

SIZE_RED = ["mlog", "sdlog", "skewlog", "kurtlog", "bowley", "tailw", "gini", "effn", "top", "theil"]
FRAC_RED = ["mean", "sd", "skew", "kurt", "med", "iqr", "bowley", "tailw"]
FEW_RED = ["mean", "sd"]
BIN_RED = ["mean"]
ORDER_RED = ["acf1", "acf2", "acf4", "acf8", "acf16", "turn"]
ETA_RED = {"B": ["eta_t", "eta_q", "eta_tq"], "J": ["eta_t"], "JQ": []}


def single_reducers(ctx, name):
    kind, cls = ctx.kind(name), ctx.cls_of(name)
    r = {"size": SIZE_RED, "frac": FRAC_RED, "few": FEW_RED, "bin": BIN_RED}[kind]
    r = list(r)
    if cls in ("B", "J"):
        r += ORDER_RED if kind in ("size", "frac") else ORDER_RED[:5]
        r += ETA_RED[cls]
    if cls == "B" and name != "tw":
        r += ["wmean"] if kind == "bin" else ["wmean", "wsd"]
    return r


def reduce_single(ctx, name, reds):
    """Values of several reducers on one vector."""
    out = {}
    kind, cls = ctx.kind(name), ctx.cls_of(name)
    y = ctx.canon(name)
    fin = np.isfinite(y)
    allfin = bool(fin.all())
    yf = y if allfin else y[fin]
    need = set(reds)
    if need & {"mlog", "sdlog", "skewlog", "kurtlog", "mean", "sd", "skew", "kurt"}:
        m, sd, sk, ku = _moments(yf)
        for a, b_ in (("mlog", m), ("sdlog", sd), ("skewlog", sk), ("kurtlog", ku),
                      ("mean", m), ("sd", sd), ("skew", sk), ("kurt", ku)):
            if a in need:
                out[a] = b_
    if need & {"med", "iqr", "bowley", "tailw", "gini"}:
        if yf.size >= 3:
            s = np.sort(yf)
            q1, q2, q3 = _q(s, .25), _q(s, .5), _q(s, .75)
            iqr = q3 - q1
            out["med"] = q2
            out["iqr"] = iqr
            out["bowley"] = (q3 + q1 - 2 * q2) / iqr if iqr > 0 else NAN
            out["tailw"] = (_q(s, .9) - _q(s, .1)) / iqr if iqr > 0 else NAN
            if "gini" in need:
                x = np.exp(s)
                k = np.arange(1, x.size + 1)
                out["gini"] = 2 * np.dot(k, x) / (x.size * x.sum()) - (x.size + 1) / x.size
        else:
            for a in ("med", "iqr", "bowley", "tailw", "gini"):
                out[a] = NAN
    if need & {"effn", "top", "theil"}:
        if yf.size >= 3:
            x = np.exp(yf)
            sx = x.sum()
            out["effn"] = sx * sx / (x.size * np.dot(x, x))
            out["top"] = x.max() / sx
            out["theil"] = np.log(sx / x.size) - yf.mean()
        else:
            out["effn"] = out["top"] = out["theil"] = NAN
    if need & {"wmean", "wsd"}:
        w = ctx.tw if allfin else ctx.tw[fin]
        if yf.size >= 3:
            wm = np.dot(w, yf) / w.sum()
            out["wmean"] = wm
            wss = np.dot(w, (yf - wm) ** 2)
            out["wsd"] = 0.0 if _zero(wss, w.sum(), wm) else np.sqrt(wss / w.sum())
        else:
            out["wmean"] = out["wsd"] = NAN
    if cls in ("B", "J"):
        g = ctx.grp(cls)
        for lag in (1, 2, 4, 8, 16):
            if "acf%d" % lag in need:
                out["acf%d" % lag] = _acf(y, g, lag)
        if "turn" in need:
            out["turn"] = _turn(y, g)
        if "eta_t" in need:
            out["eta_t"] = _eta2(yf, g if allfin else g[fin])
        if "eta_q" in need:
            out["eta_q"] = _eta2(yf, ctx.qseq if allfin else ctx.qseq[fin])
        if "eta_tq" in need:
            out["eta_tq"] = _eta2(yf, ctx.gid if allfin else ctx.gid[fin])
    return {k: out[k] for k in reds}


# ----------------------------------------------------------- pair reducers
#
#   pcor  Pearson correlation on the canonical scale
#   rcor  Pearson correlation of the ranks (ranks taken within each vector)
#   wcor  pcor weighted by the width of the blocks on the target (B only)
#   gcor  pcor after centring each vector on its mean per target sequence
#   xcor  correlation between a block's value of one vector and its
#         neighbours' value of the other, both directions pooled (B only)

def _corr_matrix(M, w=None):
    """Correlation matrix of the rows of M (all finite), optionally weighted."""
    n = M.shape[1]
    if w is None:
        mu = M.mean(axis=1, keepdims=True)
        M = M - mu
        ss = np.einsum("ij,ij->i", M, M)
        const = ss <= TOL * n * (mu[:, 0] ** 2 + 1.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            Z = M / np.sqrt(ss)[:, None]
    else:
        sw = w.sum()
        mu = (M @ w / sw)[:, None]
        M = M - mu
        ss = np.einsum("ij,ij,j->i", M, M, w)
        const = ss <= TOL * sw * (mu[:, 0] ** 2 + 1.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            Z = M * np.sqrt(w)[None, :] / np.sqrt(ss)[:, None]
    Z[const] = 0.0
    C = Z @ Z.T
    C[const, :] = np.nan
    C[:, const] = np.nan
    return C


def _center_groups(M, g):
    """Subtract from each row of M its mean within the groups g."""
    _, inv = np.unique(g, return_inverse=True)
    cnt = np.bincount(inv)
    out = np.empty_like(M)
    for k in range(M.shape[0]):
        out[k] = M[k] - (np.bincount(inv, weights=M[k]) / cnt)[inv]
    return out


def pair_corrs(arrs, w=None, g=None):
    """All pairwise correlations between 1-D arrays of equal length that may
    contain NaN; each pair uses the positions where both are finite.
    `w`: weights; `g`: groups for centring."""
    k = len(arrs)
    C = np.full((k, k), np.nan)
    if k == 0 or arrs[0].size < 3:
        return C

    def cm(rows, mask=None):
        M = np.vstack(rows) if mask is None else np.vstack([r[mask] for r in rows])
        ww = None if w is None else (w if mask is None else w[mask])
        if g is not None:
            M = _center_groups(M, g if mask is None else g[mask])
        return _corr_matrix(M, ww)

    fin = [np.isfinite(a) for a in arrs]
    full = np.array([f.all() for f in fin])
    fi = np.flatnonzero(full)
    if fi.size > 1:
        C[np.ix_(fi, fi)] = cm([arrs[j] for j in fi])
    part = np.flatnonzero(~full)
    for a in part:
        m = fin[a]
        if m.sum() < 3:
            continue
        if fi.size:
            sub = cm([arrs[a]] + [arrs[j] for j in fi], m)
            C[a, fi] = sub[0, 1:]
            C[fi, a] = sub[0, 1:]
        for b_ in part:
            if b_ > a:
                mm = m & fin[b_]
                if mm.sum() >= 3:
                    C[a, b_] = C[b_, a] = cm([arrs[a], arrs[b_]], mm)[0, 1]
    return C


def xcor_matrix(ctx, names):
    """Symmetric lag-1 cross-correlations between B vectors without NaN."""
    k = len(names)
    C = np.full((k, k), np.nan)
    i = ctx.iJ
    if i.size < 3:
        return C
    X = np.vstack([ctx.canon(v) for v in names])
    A, B = X[:, i], X[:, i + 1]
    m = (A.sum(1) + B.sum(1)) / (2 * i.size)
    A = A - m[:, None]
    B = B - m[:, None]
    ss = (np.einsum("ij,ij->i", A, A) + np.einsum("ij,ij->i", B, B)) / 2
    const = ss <= TOL * i.size * (m ** 2 + 1.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        sd = np.sqrt(ss)
        C = (A @ B.T + B @ A.T) / 2 / np.outer(sd, sd)
    C[const, :] = np.nan
    C[:, const] = np.nan
    return C


XCOR_SKIP = ("glen", "bgap", "bflip", "bsameq")


def pair_families(ctx):
    """(reducer, vector names, function returning the correlation matrix)."""
    fam = []
    B, J, JQ = ctx.names("B"), ctx.names("J"), ctx.names("JQ")
    fam.append(("pcor", "B", B, lambda: pair_corrs([ctx.canon(v) for v in B])))
    fam.append(("rcor", "B", B, lambda: pair_corrs([ctx.rank(v) for v in B])))
    fam.append(("wcor", "B", B, lambda: pair_corrs([ctx.canon(v) for v in B], w=ctx.tw)))
    fam.append(("gcor", "B", B, lambda: pair_corrs([ctx.canon(v) for v in B], g=ctx.tseq)))
    Bx = [v for v in B if v not in XCOR_SKIP]
    fam.append(("xcor", "B", Bx, lambda: xcor_matrix(ctx, Bx)))
    fam.append(("pcor", "J", J, lambda: pair_corrs([ctx.canon(v) for v in J])))
    fam.append(("rcor", "J", J, lambda: pair_corrs([ctx.rank(v) for v in J])))
    fam.append(("gcor", "J", J, lambda: pair_corrs([ctx.canon(v) for v in J], g=ctx.jseq)))
    Jq = [v for v in J if v != "sameq"]
    mix = JQ + Jq
    fam.append(("pcor", "JQ", mix, lambda: pair_corrs(
        [ctx.canon(v) for v in JQ] + [ctx.on_jq(ctx.canon(v)) for v in Jq])))
    fam.append(("rcor", "JQ", mix, lambda: pair_corrs(
        [ctx.rank(v) for v in JQ] + [ctx.on_jq(ctx.rank(v)) for v in Jq])))
    return fam


def pair_names(level, ctx):
    names = []
    for red, cls, vs, _ in pair_families(ctx):
        nfirst = len(ctx.names("JQ")) if cls == "JQ" else len(vs)
        for a in range(nfirst):
            for c in range(a + 1, len(vs)):
                names.append("%s.%s.%s~%s" % (level, red, vs[a], vs[c]))
    return names


def pair_candidates(ctx, out, level):
    for red, cls, vs, fun in pair_families(ctx):
        C = fun()
        nfirst = len(ctx.names("JQ")) if cls == "JQ" else len(vs)
        for a in range(nfirst):
            for c in range(a + 1, len(vs)):
                out["%s.%s.%s~%s" % (level, red, vs[a], vs[c])] = float(C[a, c])


def pair_one(ctx, red, a, c):
    """One pair correlation, computed alone (same value as in the matrix)."""
    ca, cc = ctx.cls_of(a), ctx.cls_of(c)
    if red == "xcor":
        return float(xcor_matrix(ctx, [a, c])[0, 1])
    getter = ctx.rank if red == "rcor" else ctx.canon
    x, y = getter(a), getter(c)
    if "JQ" in (ca, cc):
        if ca == "J":
            x = ctx.on_jq(x)
        if cc == "J":
            y = ctx.on_jq(y)
        return float(pair_corrs([x, y])[0, 1])
    w = ctx.tw if red == "wcor" else None
    g = (ctx.tseq if ca == "B" else ctx.jseq) if red == "gcor" else None
    return float(pair_corrs([x, y], w=w, g=g)[0, 1])


# ------------------------------------------------ matrix and permutations

MATRIX = ["nmi", "hcond", "cramer", "effpart", "second", "sv1", "cca1"]
PERM = ["foot", "fix", "desc", "cyc", "lis", "rho", "rec"]


def matrix_candidates(ctx, want=None):
    """Statistics of the table of aligned target bases by pair of sequences."""
    out = {}
    tu, ti = np.unique(ctx.tseq, return_inverse=True)
    qu, qi = np.unique(ctx.qseq, return_inverse=True)
    W = np.zeros((tu.size, qu.size))
    np.add.at(W, (ti, qi), ctx.tw)
    tot = W.sum()
    if tot <= 0:
        return {k: NAN for k in MATRIX}
    Pm = W / tot
    pr, pc = Pm.sum(1), Pm.sum(0)
    with np.errstate(divide="ignore", invalid="ignore"):
        H = lambda p: -np.sum(np.where(p > 0, p * np.log(p), 0.0))
        hr, hc, hj = H(pr), H(pc), H(Pm)
        mi = hr + hc - hj
        out["nmi"] = mi / np.sqrt(hr * hc) if hr > 0 and hc > 0 else NAN
        out["hcond"] = (hj - hr) / np.log(qu.size) if qu.size > 1 else 0.0
        E = np.outer(pr, pc)
        chi = np.sum(np.where(E > 0, (Pm - E) ** 2 / E, 0.0))
        k = min(tu.size, qu.size)
        out["cramer"] = np.sqrt(chi / (k - 1)) if k > 1 else NAN
        R = Pm / pr[:, None]
        hrow = -np.sum(np.where(R > 0, R * np.log(R), 0.0), axis=1)
        out["effpart"] = float(np.dot(pr, np.exp(hrow)))
        if qu.size > 1:
            out["second"] = float(np.dot(pr, np.sort(R, axis=1)[:, -2]))
        else:
            out["second"] = 0.0
        sv = np.linalg.svd(Pm, compute_uv=False)
        out["sv1"] = float(sv[0] ** 2 / np.sum(sv ** 2))
        S = Pm / np.sqrt(E)
        S[~np.isfinite(S)] = 0
        sv = np.linalg.svd(S, compute_uv=False)
        out["cca1"] = float(sv[1]) if sv.size > 1 else NAN
    return out


def _lis_len(seq):
    tails = []
    for x in seq:
        k = bisect_left(tails, x)
        if k == len(tails):
            tails.append(x)
        else:
            tails[k] = x
    return len(tails)


def perm_candidates(ctx, want=None):
    """Statistics of the permutation that maps the target order of the blocks
    to their query order, within each pair of sequences, read in the
    orientation of the majority strand.  Averages are weighted by the number
    of blocks; pairs of sequences with a single block are left out."""
    out = {k: NAN for k in PERM}
    n = ctx.n
    if n < 3:
        return out
    gid = ctx.gid
    o_t = np.lexsort((np.arange(n), gid))            # target order within group
    g = gid[o_t]
    cnt = np.bincount(gid)
    start = np.r_[0, np.cumsum(cnt)[:-1]]
    k = np.arange(n) - start[g]                      # target rank within group
    o_q = np.lexsort((ctx.qrank, gid))
    qr = np.empty(n, dtype=np.int64)
    qr[o_q] = np.arange(n) - start[gid[o_q]]         # query rank within group, by block
    pi = qr[o_t]
    ng = cnt[g]
    pi = np.where(ctx.gmaj[g] > 0, pi, ng - 1 - pi)
    big = ng >= 2
    if not big.any():
        return out
    wsum = big.sum()
    d = np.abs(pi - k).astype(float)
    foot_g = np.bincount(g, weights=d, minlength=cnt.size) / np.maximum(np.floor(cnt * cnt / 2.0), 1)
    out["foot"] = float(np.dot(foot_g[cnt >= 2], cnt[cnt >= 2]) / wsum)
    out["fix"] = float(((pi == k) & big).sum() / wsum)
    same = g[1:] == g[:-1]
    out["desc"] = float(((pi[1:] < pi[:-1]) & same).sum() / max(same.sum(), 1))
    d2 = np.bincount(g, weights=d * d, minlength=cnt.size)
    c = cnt.astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rho_g = 1 - 6 * d2 / (c * (c * c - 1))
    out["rho"] = float(np.dot(rho_g[cnt >= 2], cnt[cnt >= 2]) / wsum)
    P = start[g] + pi                                # global permutation
    # records: left-to-right maxima and right-to-left minima, averaged so that
    # the count does not depend on the orientation of the sequences
    rec = (P == np.maximum.accumulate(P)).astype(float) \
        + (P == np.minimum.accumulate(P[::-1])[::-1])
    out["rec"] = float((rec * big).sum() / 2 / wsum)
    if want is None or "cyc" in want:
        A = coo_matrix((np.ones(n, dtype=np.int8), (np.arange(n), P)), shape=(n, n))
        ncomp, lab = connected_components(A, directed=True, connection="weak")
        first = np.unique(lab, return_index=True)[1]
        cyc_g = np.bincount(g[first], minlength=cnt.size)
        out["cyc"] = float(cyc_g[cnt >= 2].sum() / wsum)
    if want is None or "lis" in want:
        tot = 0
        pl = pi.tolist()
        st = start.tolist()
        for gi in np.flatnonzero(cnt >= 2):
            s = pl[st[gi]: st[gi] + int(cnt[gi])]
            tot += _lis_len(s)
        out["lis"] = tot / wsum
    return out


# ------------------------------------------------------- ordinal patterns

PAT = ["o3A", "o3B", "o3C", "o3D", "H3", "H4", "s3aaa", "s3bbb", "s3aba", "s3bab",
       "s3aab", "s3abb", "Hs3", "Hjoint"]


def _entropy(counts, k):
    tot = counts.sum()
    if tot <= 0:
        return NAN
    p = counts[counts > 0] / tot
    return float(-np.sum(p * np.log(p)) / np.log(k))


_M3 = np.array([0, 2, 1, 3, 4, 6, 5, 7])                      # mirror of an ordinal code of 3
_S3 = np.array([4 * (c & 1) + (c & 2) + (c >> 2) for c in range(8)])   # mirror of a strand code
_M4 = np.array([sum((((c >> src) & 1) << dst) for src, dst in
                    ((0, 5), (5, 0), (1, 4), (4, 1), (2, 2), (3, 3))) for c in range(64)])


def pat_candidates(ctx):
    """Frequencies and entropies of local patterns of consecutive blocks.

    Ordinal patterns: relative order on the query of 3 (or 4) consecutive
    blocks that are on the same target and the same query sequence, read in
    the orientation of the majority strand.
    Strand patterns: strands of 3 consecutive blocks on the same target
    sequence, a = majority strand of the block's pair of sequences, b = the
    other.
    A pattern and its mirror image (the same blocks read from the other end
    of the sequence) are counted as one class, because the reading direction
    is arbitrary: A = 123, B = 321, C = 132 or 213, D = 231 or 312; aab or
    baa; abb or bba.  Entropies are divided by the logarithm of the number
    of classes (4, 16, 6 and 28)."""
    out = {k: NAN for k in PAT}
    n = ctx.n
    if n < 4:
        return out
    t, g = ctx.tseq, ctx.gid
    r = ctx.qrank * ctx.gmaj[g]                       # oriented query rank
    sm = (ctx.vec("smaj") < 0).astype(np.int64)        # 0 = majority, 1 = minority
    s3 = 4 * sm[:-2] + 2 * sm[1:-1] + sm[2:]
    ok3 = (g[2:] == g[:-2]) & (g[1:-1] == g[:-2])
    if ok3.sum() >= 3:
        r0, r1, r2 = r[:-2][ok3], r[1:-1][ok3], r[2:][ok3]
        code = ((r0 < r1) + 2 * (r1 < r2) + 4 * (r0 < r2)).astype(np.int64)
        c = np.bincount(np.minimum(code, _M3[code]), minlength=8).astype(float)
        tot = c.sum()
        # 123 -> 7, 321 -> 0, 132 or 213 -> 5, 231 or 312 -> 1
        out["o3A"], out["o3B"], out["o3C"], out["o3D"] = c[7] / tot, c[0] / tot, c[5] / tot, c[1] / tot
        out["H3"] = _entropy(c, 4)
        sc = s3[ok3]
        joint = 8 * code + sc
        mirror = 8 * _M3[code] + _S3[sc]
        out["Hjoint"] = _entropy(np.bincount(np.minimum(joint, mirror), minlength=64).astype(float), 28)
    ok4 = (g[3:] == g[:-3]) & (g[2:-1] == g[:-3]) & (g[1:-2] == g[:-3])
    if ok4.sum() >= 3:
        q = [r[k: n - 3 + k][ok4] for k in range(4)]
        code = np.zeros(q[0].size, dtype=np.int64)
        bit = 0
        for a in range(4):
            for b_ in range(a + 1, 4):
                code += (q[a] < q[b_]).astype(np.int64) << bit
                bit += 1
        out["H4"] = _entropy(np.bincount(np.minimum(code, _M4[code]), minlength=64).astype(float), 16)
    oks = t[2:] == t[:-2]
    if oks.sum() >= 3:
        code = s3[oks]
        c = np.bincount(np.minimum(code, _S3[code]), minlength=8).astype(float)
        tot = c.sum()
        # aaa -> 0, bbb -> 7, aba -> 2, bab -> 5, aab or baa -> 1, abb or bba -> 3
        out["s3aaa"], out["s3bbb"], out["s3aba"], out["s3bab"] = c[0] / tot, c[7] / tot, c[2] / tot, c[5] / tot
        out["s3aab"], out["s3abb"] = c[1] / tot, c[3] / tot
        out["Hs3"] = _entropy(c, 6)
    return out


# ---------------------------------------------------------------- features

def _empty_ctx(level):
    b = dict(ts=np.zeros(0), te=np.zeros(0), qs=np.zeros(0), qe=np.zeros(0),
             tseq=np.zeros(0, int), qseq=np.zeros(0, int), strand=np.zeros(0),
             tL=np.zeros(1), qL=np.zeros(1))
    return Ctx(b, level)


def candidate_names(level):
    """All candidate names of one level (no data needed)."""
    ctx = _empty_ctx(level)
    names = []
    for cls in ("B", "J", "JQ"):
        for v in ctx.names(cls):
            for r in single_reducers(ctx, v):
                names.append("%s.%s.%s" % (level, r, v))
    names += pair_names(level, ctx)
    names += ["%s.mat.%s" % (level, m) for m in MATRIX]
    names += ["%s.perm.%s" % (level, m) for m in PERM]
    names += ["%s.pat.%s" % (level, m) for m in PAT]
    return names


PAIR_RED = ("pcor", "rcor", "wcor", "gcor", "xcor")


def features(b, level, want=None):
    """All candidate values for one alignment at one level, or only those in
    `want` (a collection of candidate names)."""
    ctx = Ctx(b, level)
    out = {}
    if want is None:
        for cls in ("B", "J", "JQ"):
            for v in ctx.names(cls):
                reds = single_reducers(ctx, v)
                res = reduce_single(ctx, v, reds)
                for r in reds:
                    out["%s.%s.%s" % (level, r, v)] = res[r]
        pair_candidates(ctx, out, level)
        for k, v in matrix_candidates(ctx).items():
            out["%s.mat.%s" % (level, k)] = v
        for k, v in perm_candidates(ctx).items():
            out["%s.perm.%s" % (level, k)] = v
        for k, v in pat_candidates(ctx).items():
            out["%s.pat.%s" % (level, k)] = v
        return out
    singles = {}
    mat = perm = pat = None
    for name in want:
        lev, red, rest = name.split(".", 2)
        if lev != level:
            continue
        if red == "mat":
            if mat is None:
                mat = matrix_candidates(ctx)
            out[name] = mat[rest]
        elif red == "perm":
            if perm is None:
                perm = perm_candidates(ctx, want={w.split(".", 2)[2] for w in want
                                                   if w.split(".", 2)[1] == "perm"})
            out[name] = perm[rest]
        elif red == "pat":
            if pat is None:
                pat = pat_candidates(ctx)
            out[name] = pat[rest]
        elif red in PAIR_RED:
            a, c = rest.split("~")
            out[name] = pair_one(ctx, red, a, c)
        else:
            singles.setdefault(rest, []).append(red)
    for v, reds in singles.items():
        res = reduce_single(ctx, v, reds)
        for r in reds:
            out["%s.%s.%s" % (level, r, v)] = res[r]
    return out


# ------------------------------------------------- arrangement-free candidates

INTRINSIC = {"tw", "qw", "asym", "aasym", "ident", "mismf", "gapf", "sdens", "gcol",
             "alen", "glen", "nmem", "tdens", "qdens"}
MARGINAL_RED = set(SIZE_RED) | set(FRAC_RED)


def arrangement_free(name):
    """TRUE for candidates that, by construction, do not use the position,
    the order or the strand of the blocks: descriptive statistics of block
    attributes and correlations between block attributes.  Their value is
    the same in any arrangement of the same blocks."""
    lev, red, rest = name.split(".", 2)
    if red in ("pcor", "rcor", "wcor"):
        a, c = rest.split("~")
        return a in INTRINSIC and c in INTRINSIC
    if red in MARGINAL_RED:
        return rest in INTRINSIC or rest == "tgap"
    if red in ("wmean", "wsd"):
        return rest in INTRINSIC
    return False
