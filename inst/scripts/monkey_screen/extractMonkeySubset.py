#!/usr/bin/env python3
"""
extractMonkeySubset.py - select and export a subset of aligned genome pairs.

For each Branch found in both AlignedGenomes/<version>/<seed>/ and
StatsOnPairs/latest/:

 1. read the .o2o.matrix.txt file of every ordered pair and compute the
    P distance (1 - identical / aligned ACGT pairs; "L3" in GenomicBreaks);
 2. build a UPGMA tree on the symmetrised distances;
 3. cut the tree into K clades and pick genomes in each clade so that the
    selection contains close, intermediate and distant pairs;
 4. check that every ordered pair of the selection has a blasttabplus file,
    a matrix file, and a YAML statistics file that is not older than the
    alignment and contains the keys of the current statistics pipeline;
 5. report sizes (default), or with --write make one tar file per Branch.

Standard library only (Python >= 3.6).  Paths are built from the directory
names; no recursive listing is done.  Reads many small files on /bucket: run
it on a compute node (/bucket is read-only there, which is enough) and write
the output to /flash.

Example:
  srun -p compute -c 1 --mem 4G -t 0-4 --pty \\
    python3 extractMonkeySubset.py --out /flash/LuscombeU/$USER/monkey
  # read summary.tsv, choose Branches, then:
  srun -p compute -c 1 --mem 4G -t 0-4 --pty \\
    python3 extractMonkeySubset.py --out /flash/LuscombeU/$USER/monkey \\
      --write --branches Primates,Buchnera
"""
import argparse
import os
import re
import sys
import tarfile

ROOT = "/bucket/LuscombeU/common/ScramblingInTheTree"
ACGT = "ACGT"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", default=ROOT)
    p.add_argument("--version", default="3.0", help="AlignedGenomes version directory")
    p.add_argument("--seeds", default="RY4,MAM8")
    p.add_argument("--stats", default="latest", help="StatsOnPairs subdirectory")
    p.add_argument("--branches", default="", help="comma-separated Branch names (default: all eligible)")
    p.add_argument("--out", required=True, help="output directory (on /flash)")
    p.add_argument("--clades", type=int, default=4, help="number of clades cut from the tree")
    p.add_argument("--genomes", type=int, default=12, help="genomes selected per Branch")
    p.add_argument("--min-dist", type=float, default=0.002,
                   help="smallest P distance allowed between two selected genomes")
    p.add_argument("--yaml-keys", default="index_strandRand_target_5,P",
                   help="keys that a YAML file must contain to count as up to date")
    p.add_argument("--allow-stale", action="store_true",
                   help="do not require YAML files to be newer than the alignments")
    p.add_argument("--threads", type=int, default=8, help="parallel reads of matrix files")
    p.add_argument("--write", action="store_true", help="write the tar files (default: report only)")
    return p.parse_args()


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- distances

def p_distance(path):
    """P distance from a .o2o.matrix.txt file; None if missing or empty."""
    try:
        with open(path) as fh:
            rows = [l.split() for l in fh if l.strip() and not l.startswith("#")]
    except OSError:
        return None
    if len(rows) < 5:
        return None
    cols = rows[0]
    if len(cols) < 4:
        return None
    same = total = 0.0
    try:
        for r in rows[1:]:
            if r[0] not in ACGT or len(r) != len(cols) + 1:
                continue
            for letter, v in zip(cols, r[1:]):
                if letter in ACGT:
                    total += float(v)
                    if letter == r[0]:
                        same += float(v)
    except ValueError:
        return None
    if total == 0:
        return None
    return 1.0 - same / total


def branch_distances(aln_dir, ids, cache, threads):
    """Ordered-pair P distances, cached in a 3-column file.  The cache is
    appended to after each target, so an interrupted run resumes from it."""
    from concurrent.futures import ThreadPoolExecutor
    d = {}
    if os.path.exists(cache):
        with open(cache) as fh:
            for l in fh:
                f = l.rstrip("\n").split("\t")
                if len(f) == 3 and l.endswith("\n"):       # skip a truncated last line
                    d[(f[0], f[1])] = None if f[2] == "NA" else float(f[2])
        with open(cache, "w") as fh:                 # rewrite without damaged lines
            for (t, q), v in d.items():
                fh.write("%s\t%s\t%s\n" % (t, q, "NA" if v is None else repr(v)))
    path = lambda t, q: os.path.join(
        aln_dir, "results_" + t, "alignment", t + "___" + q + ".o2o.matrix.txt")
    with ThreadPoolExecutor(threads) as pool, open(cache, "a") as fh:
        for n, t in enumerate(ids, 1):
            todo = [q for q in ids if q != t and (t, q) not in d]
            if not todo:
                continue
            for q, v in zip(todo, pool.map(lambda q: p_distance(path(t, q)), todo)):
                d[(t, q)] = v
                fh.write("%s\t%s\t%s\n" % (t, q, "NA" if v is None else repr(v)))
            fh.flush()
            log("    matrices read for target %d/%d" % (n, len(ids)))
    return d


def symmetrise(d, ids):
    """Mean of both orders; one order if the other is missing; None if both are."""
    s = {}
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            v = [x for x in (d.get((a, b)), d.get((b, a))) if x is not None]
            s[(a, b)] = s[(b, a)] = sum(v) / len(v) if v else None
    return s


# --------------------------------------------------------------------- tree

def upgma(ids, s):
    """Average-linkage clustering.  Returns the list of merges
    (members_left, members_right, height) in order of increasing height."""
    clusters = {i: [g] for i, g in enumerate(ids)}
    dist = {}
    keys = list(clusters)
    for x, i in enumerate(keys):
        for j in keys[x + 1:]:
            dist[(i, j)] = s[(ids[i], ids[j])]
    merges, nxt = [], len(ids)
    while len(clusters) > 1:
        (i, j), h = min(dist.items(), key=lambda kv: (kv[1], kv[0]))
        ni, nj = len(clusters[i]), len(clusters[j])
        merges.append((clusters[i], clusters[j], h))
        new = {}
        for k in clusters:
            if k in (i, j):
                continue
            dik = dist[(min(i, k), max(i, k))]
            djk = dist[(min(j, k), max(j, k))]
            new[k] = (dik * ni + djk * nj) / (ni + nj)
        for k in list(dist):
            if i in k or j in k:
                del dist[k]
        members = clusters.pop(i) + clusters.pop(j)
        for k, v in new.items():
            dist[(k, nxt)] = v
        clusters[nxt] = members
        nxt += 1
    return merges


def newick(ids, merges):
    node = {g: (g, 0.0) for g in ids}        # representative member -> (text, height)
    rep = {g: g for g in ids}                # member -> current cluster key
    for left, right, h in merges:
        a, b = rep[left[0]], rep[right[0]]
        (ta, ha), (tb, hb) = node.pop(a), node.pop(b)
        txt = "(%s:%.6f,%s:%.6f)" % (ta, max(h / 2 - ha, 0), tb, max(h / 2 - hb, 0))
        node[a] = (txt, h / 2)
        for g in left + right:
            rep[g] = a
    return list(node.values())[0][0] + ";"


def cut(ids, merges, k):
    """Clades obtained by undoing the last k-1 merges."""
    k = max(1, min(k, len(ids)))
    clades = [list(ids)] if not merges else [merges[-1][0] + merges[-1][1]]
    for left, right, _ in reversed(merges[len(merges) - (k - 1):] if k > 1 else []):
        whole = set(left) | set(right)
        clades = [c for c in clades if set(c) != whole] + [left, right]
    return clades


# ---------------------------------------------------------------- selection

def pick_in_clade(clade, n, s, min_dist, chosen):
    """Seed = member closest on average to its clade; then its nearest
    neighbour beyond min_dist; then farthest-point additions."""
    ok = lambda g, sel: all(s[(g, x)] >= min_dist for x in sel)
    pool = [g for g in sorted(clade) if ok(g, chosen)]
    sel = []
    if not pool or n < 1:
        return sel
    mean = lambda g: sum(s[(g, x)] for x in clade if x != g) / max(len(clade) - 1, 1)
    sel.append(min(pool, key=lambda g: (mean(g), g)))
    while len(sel) < n:
        cand = [g for g in pool if g not in sel and ok(g, sel)]
        if not cand:
            break
        if len(sel) == 1:
            sel.append(min(cand, key=lambda g: (s[(g, sel[0])], g)))
        else:
            sel.append(max(cand, key=lambda g: (min(s[(g, x)] for x in sel), g)))
    return sel


def select(ids, s, n_clades, n_genomes, min_dist):
    merges = upgma(ids, s)
    clades = sorted(cut(ids, merges, n_clades), key=lambda c: (-len(c), sorted(c)[0]))
    quota = [n_genomes // len(clades)] * len(clades)
    for i in range(n_genomes - sum(quota)):
        quota[i] += 1
    chosen, origin = [], {}
    for _ in range(3):                       # redistribute what small clades cannot give
        chosen, origin = [], {}
        for ci, (clade, q) in enumerate(zip(clades, quota)):
            got = pick_in_clade(clade, q, s, min_dist, chosen)
            chosen += got
            origin.update((g, ci + 1) for g in got)
            quota[ci] = len(got)
        deficit = n_genomes - len(chosen)
        if deficit <= 0:
            break
        for ci in sorted(range(len(clades)), key=lambda i: quota[i] - len(clades[i])):
            if deficit and quota[ci] < len(clades[ci]):
                quota[ci] += 1
                deficit -= 1
        if deficit == n_genomes - len(chosen):
            break
    return chosen, origin, merges


# ------------------------------------------------------------------- checks

def removed_patterns(*dirs):
    pats = []
    for d in dirs:
        f = os.path.join(d, "removedAssemblies.txt")
        if os.path.exists(f):
            with open(f) as fh:
                for l in fh:
                    pat = l.split("\t")[0].strip()
                    if pat:
                        try:
                            pats.append(re.compile(pat))
                        except re.error:
                            pats.append(re.compile(re.escape(pat)))
    return pats


def pair_files(aln_dir, stat_dir, t, q):
    base = os.path.join(aln_dir, "results_" + t, "alignment", t + "___" + q)
    return {"blasttabplus": base + ".blasttabplus.gz",
            "matrix": base + ".o2o.matrix.txt",
            "o2o": base + ".o2o.tsv",
            "yaml": os.path.join(stat_dir, "genomicbreaks", t + "___" + q + ".yaml")}


def check_pair(files, keys, allow_stale):
    """Returns (problem or None, size in bytes)."""
    size = 0
    for kind in ("blasttabplus", "matrix", "yaml"):
        try:
            st = os.stat(files[kind])
        except OSError:
            return "missing " + kind, 0
        if st.st_size == 0:
            return "empty " + kind, 0
        size += st.st_size
    if not allow_stale and os.stat(files["yaml"]).st_mtime < os.stat(files["blasttabplus"]).st_mtime:
        return "yaml older than alignment", size
    with open(files["yaml"], errors="replace") as fh:
        txt = fh.read()
    for k in keys:
        if not re.search(r"(^|[^A-Za-z0-9_])" + re.escape(k) + r"\s*:", txt):
            return "yaml lacks key " + k, size
    return None, size


# --------------------------------------------------------------------- main

def process(branch, seed, a, keys):
    aln_dir = os.path.join(a.root, "AlignedGenomes", a.version, seed, branch)
    stat_dir = os.path.join(a.root, "StatsOnPairs", a.stats, branch)
    out = os.path.join(a.out, seed + "_" + branch)
    os.makedirs(out, exist_ok=True)
    ids = sorted(e[len("results_"):] for e in os.listdir(aln_dir) if e.startswith("results_"))
    pats = removed_patterns(aln_dir, stat_dir)
    removed = [g for g in ids if any(p.search(g) for p in pats)]
    ids = [g for g in ids if g not in removed]
    log("  %d genomes (%d excluded by removedAssemblies.txt)" % (len(ids), len(removed)))
    if len(ids) < 4:
        return dict(status="too few genomes")
    d = branch_distances(aln_dir, ids, os.path.join(out, "P_ordered_pairs.tsv"), a.threads)
    s = symmetrise(d, ids)
    # drop genomes with undefined distances, worst first
    while True:
        bad = {g: sum(1 for x in ids if x != g and s[(g, x)] is None) for g in ids}
        worst = max(ids, key=lambda g: (bad[g], g))
        if bad[worst] == 0:
            break
        ids.remove(worst)
        removed.append(worst)
    excluded, problems = set(), []
    while True:
        cand = [g for g in ids if g not in excluded]
        if len(cand) < 4:
            return dict(status="no complete subset", detail="; ".join(problems[:3]))
        chosen, origin, merges = select(cand, s, a.clades, a.genomes, a.min_dist)
        fail, total, stale = {}, 0, None
        for t in chosen:
            for q in chosen:
                if t == q:
                    continue
                why, size = check_pair(pair_files(aln_dir, stat_dir, t, q), keys, a.allow_stale)
                total += size
                if why:
                    if why.startswith("yaml"):
                        stale = why
                    fail[t] = fail.get(t, 0) + 1
                    fail[q] = fail.get(q, 0) + 1
                    problems.append("%s___%s: %s" % (t, q, why))
        if not fail:
            break
        # a Branch whose statistics are out of date fails everywhere: stop early
        if stale and len(fail) == len(chosen) and min(fail.values()) >= len(chosen) - 1:
            return dict(status="statistics not up to date", detail=stale)
        excluded.add(max(fail, key=lambda g: (fail[g], g)))
    with open(os.path.join(out, "tree_all_genomes.nwk"), "w") as fh:
        fh.write(newick(cand, merges) + "\n")
    with open(os.path.join(out, "selection.tsv"), "w") as fh:
        fh.write("id\tclade\n")
        for g in chosen:
            fh.write("%s\t%d\n" % (g, origin[g]))
    with open(os.path.join(out, "excluded.txt"), "w") as fh:
        for g in removed:
            fh.write(g + "\tremovedAssemblies.txt or undefined distances\n")
        for g in sorted(excluded):
            fh.write(g + "\tincomplete files in a selection\n")
        for p in problems:
            fh.write("# " + p + "\n")
    dists = sorted(s[(x, y)] for i, x in enumerate(chosen) for y in chosen[i + 1:])
    res = dict(status="ok", genomes=len(chosen), pairs=len(chosen) * (len(chosen) - 1),
               bytes=total, pmin=dists[0], pmed=dists[len(dists) // 2], pmax=dists[-1])
    if a.write:
        tar_path = os.path.join(a.out, seed + "_" + branch + ".tar")
        with tarfile.open(tar_path, "w") as tar:
            top = seed + "_" + branch
            for f in ("P_ordered_pairs.tsv", "tree_all_genomes.nwk", "selection.tsv", "excluded.txt"):
                tar.add(os.path.join(out, f), arcname=top + "/" + f)
            for f in (branch + ".tsv", "removedAssemblies.txt"):
                if os.path.exists(os.path.join(aln_dir, f)):
                    tar.add(os.path.join(aln_dir, f), arcname=top + "/" + f)
            for t in chosen:
                for q in chosen:
                    if t != q:
                        for kind, f in pair_files(aln_dir, stat_dir, t, q).items():
                            if os.path.exists(f):
                                tar.add(f, arcname=top + "/" + kind + "/" + os.path.basename(f))
        res["tar"] = tar_path
    return res


def main():
    a = parse_args()
    keys = [k for k in a.yaml_keys.split(",") if k]
    wanted = set(b for b in a.branches.split(",") if b)
    os.makedirs(a.out, exist_ok=True)
    stats_root = os.path.join(a.root, "StatsOnPairs", a.stats)
    with_stats = set(os.listdir(stats_root))
    rows = []
    for seed in a.seeds.split(","):
        seed_dir = os.path.join(a.root, "AlignedGenomes", a.version, seed)
        if not os.path.isdir(seed_dir):
            continue
        for branch in sorted(os.listdir(seed_dir)):
            if wanted and branch not in wanted:
                continue
            if not os.path.isdir(os.path.join(seed_dir, branch)):
                continue
            log("%s/%s" % (seed, branch))
            if branch not in with_stats:
                res = dict(status="no StatsOnPairs directory")
            else:
                try:
                    res = process(branch, seed, a, keys)
                except Exception as e:           # keep going with the other Branches
                    res = dict(status="error", detail="%s: %s" % (type(e).__name__, e))
            log("  -> " + res["status"] + (" (" + res["detail"] + ")" if res.get("detail") else ""))
            rows.append((seed, branch, res))
    with open(os.path.join(a.out, "summary.tsv"), "w") as fh:
        fh.write("seed\tbranch\tstatus\tgenomes\tpairs\tMB\tP_min\tP_median\tP_max\tdetail\n")
        for seed, branch, r in rows:
            fh.write("\t".join([seed, branch, r["status"], str(r.get("genomes", "")),
                                str(r.get("pairs", "")),
                                "%.1f" % (r["bytes"] / 1e6) if "bytes" in r else "",
                                "%.4f" % r["pmin"] if "pmin" in r else "",
                                "%.4f" % r["pmed"] if "pmed" in r else "",
                                "%.4f" % r["pmax"] if "pmax" in r else "",
                                r.get("detail", "")]) + "\n")
    log("Summary: " + os.path.join(a.out, "summary.tsv"))
    log("Total: %.1f MB in %d Branches" % (
        sum(r.get("bytes", 0) for _, _, r in rows) / 1e6,
        sum(1 for _, _, r in rows if r["status"] == "ok")))


if __name__ == "__main__":
    main()
