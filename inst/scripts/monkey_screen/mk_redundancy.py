"""Transferable multivariate redundancy.

Can a candidate be predicted from all the reference statistics together?
Candidate and references are first turned into ranks within each Branch
(so that differences between Branches, which identify the genomes, carry no
information).  A model is fitted on all Branches but one and predicts the
held-out Branch; the score is the median, over held-out Branches, of the
squared rank correlation between prediction and truth (0 when negative).

ridge():  linear in the ranks, all candidates at once, best of several
          penalties per candidate (favours redundancy: conservative).
trees():  gradient-boosted trees (non-linear, interactions), one candidate.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.stats import rankdata, spearmanr


def within_rank(X, branch):
    """Ranks within Branch scaled to (0, 1); NaN stays NaN."""
    R = np.full(X.shape, np.nan)
    for b in np.unique(branch):
        m = np.flatnonzero(branch == b)
        for j in range(X.shape[1]):
            x = X[m, j]
            f = np.isfinite(x)
            if f.sum() >= 2:
                R[m[f], j] = rankdata(x[f]) / (f.sum() + 1)
    return R


def _score(pred, Y, branch):
    """Median over Branches of squared positive rank correlation, per column."""
    out = []
    for b in np.unique(branch):
        m = branch == b
        r = np.zeros(Y.shape[1])
        for j in range(Y.shape[1]):
            y, p = Y[m, j], pred[m, j]
            f = np.isfinite(y)
            if f.sum() >= 12 and np.std(y[f]) > 0 and np.std(p[f]) > 0:
                c = spearmanr(p[f], y[f])[0]
                r[j] = c * c if c > 0 else 0.0
            else:
                r[j] = np.nan
        out.append(r)
    with np.errstate(all="ignore"):
        return np.nanmedian(np.array(out), axis=0)


def ridge(Y, X, branch, lambdas=(0.01, 0.1, 1.0, 10.0)):
    """Y: candidates (rows x C), X: references (rows x R), both as within-Branch
    ranks.  Returns the best score per candidate over the penalties."""
    Xf = np.where(np.isfinite(X), X, 0.5) - 0.5
    Yf = np.where(np.isfinite(Y), Y, 0.5) - 0.5
    best = np.zeros(Y.shape[1])
    for lam in lambdas:
        pred = np.zeros(Y.shape)
        for b in np.unique(branch):
            te = branch == b
            A = Xf[~te]
            G = A.T @ A + lam * A.shape[0] / 100.0 * np.eye(A.shape[1])
            beta = np.linalg.solve(G, A.T @ Yf[~te])
            pred[te] = Xf[te] @ beta
        best = np.fmax(best, _score(pred, Y, branch))
    return best


def trees(y, X, branch, seed=0):
    """y: one candidate as within-Branch ranks (NaN allowed)."""
    from sklearn.ensemble import HistGradientBoostingRegressor
    pred = np.full(y.size, np.nan)
    f = np.isfinite(y)
    for b in np.unique(branch):
        te = branch == b
        m = HistGradientBoostingRegressor(max_iter=100, learning_rate=0.1, max_leaf_nodes=15,
                                          min_samples_leaf=20, random_state=seed)
        m.fit(X[~te & f], y[~te & f])
        pred[te] = m.predict(X[te])
    return float(_score(pred[:, None], y[:, None], branch)[0])
