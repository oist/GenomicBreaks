import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import mk_common as mk, mk_engine as en, mk_null as nu

rng = np.random.default_rng(1)
for br, k in (("MAM8_Saccharomycotina", 7), ("RY4_Apoidea", 5), ("RY4_Eupercaria", 40), ("MAM8_Buchnera", 3)):
    t, q = mk.pairs_of(br)[k]
    b = mk.load_cache(br, t, q)
    b.pop("tnames"); b.pop("qnames")
    c = mk.coalesce(b)
    t0 = time.time(); fr = en.features(b, "raw"); t1 = time.time(); fc = en.features(c, "coa"); t2 = time.time()
    names = en.candidate_names("raw") + en.candidate_names("coa")
    f = {**fr, **fc}
    assert set(names) == set(f), (set(names) ^ set(f))
    vals = np.array([f[x] for x in names], float)
    print(br, "blocks", len(b["ts"]), "chains", len(c["ts"]), "candidates", len(names),
          "NaN", int(np.isnan(vals).sum()), "time raw %.1fs coa %.1fs" % (t1 - t0, t2 - t1))
    # invariance to representation
    b2 = nu.rerepresent(b, rng); c2 = mk.coalesce(b2)
    g = {**en.features(b2, "raw"), **en.features(c2, "coa")}
    v2 = np.array([g[x] for x in names], float)
    both = np.isfinite(vals) & np.isfinite(v2)
    diff = np.abs(vals - v2) / np.maximum(np.abs(vals), 1e-9)
    bad = [(names[i], vals[i], v2[i]) for i in np.flatnonzero((both & (diff > 1e-6)) | (np.isfinite(vals) != np.isfinite(v2)))]
    from collections import Counter
    print("  not invariant:", len(bad), Counter(x[0].split(".")[1] for x in bad).most_common(8))
    print("     ", [(x[0], round(float(x[1]),6), round(float(x[2]),6)) for x in bad[:10]])
    # consistency of the `want` path
    sub = list(rng.choice(names, 200, replace=False))
    h = {**en.features(b, "raw", want=sub), **en.features(c, "coa", want=sub)}
    d = max(abs(h[x] - f[x]) for x in sub if np.isfinite(f[x]))
    print("  want-path max abs difference:", d, "missing", [x for x in sub if x not in h][:3])
    # swap symmetry of coalescing
    sw = mk.swap(b); cs = mk.coalesce(sw)
    print("  chains after swap:", len(cs["ts"]))
    # shuffle
    sh = nu.shuffle(b, rng)
    print("  shuffle: chains", len(mk.coalesce(sh)["ts"]), "qw kept", bool(((sh["qe"]-sh["qs"]) == (b["qe"]-b["qs"])).all()))
