# AGENTS.md — Orientation for GenomicBreaks

This file orients contributors and AI coding sessions working on
`GenomicBreaks`. For a detailed review of known engineering gaps, see
[`CODE_REVIEW.md`](CODE_REVIEW.md).

## Status / urgency

**There is no official public release of this package yet.** Consequently, the
engineering issues recorded in `docs/code_review.md` (license placeholder,
missing CI, unpinned dependencies, etc.) are **low urgency**: fix them
opportunistically, not as blockers. Do not let them stall feature work.

The detailed reusability/code-review notes are in `CODE_REVIEW.md` (repo root),
not `docs/code_review.md`.

## Commit Messages

When creating commits, include the line `AI-generated: yes` in the commit
message body (not the subject line).

## What this package is

An R package built on the **Bioconductor** ecosystem (but not part of it) for
analysing **pairwise whole-genome alignments** where evolution has scrambled
gene order. It loads alignments (GFF3 / MAF), coalesces colinear fragments,
flags structural variants (inversions, translocations, double/twin inversions),
computes synteny and nucleotide-distance metrics, and plots macrosynteny /
Oxford / pair-of-chromosome figures.

## Core data model

- `GBreaks` is an S4 class wrapping `GenomicRanges::GRanges`
  (`R/GenomicBreaks.R:22`).
- Main ranges = **target** genome coordinates.
- `query` metadata column (a `GRanges`) = **query** genome coordinates.
- Strand is carried on the target ranges; query ranges stay on the plus strand
  by convention.

## Repository layout

| Path | Contents |
|------|----------|
| `R/` | R source, one topic per file (see `Collate:` in `DESCRIPTION`). |
| `src/` | Rcpp C++ MAF reader (`readMaf.cpp`) + `RcppExports.cpp`; partly from LAST (GPLv3+). |
| `data/` | Lazy-loaded example datasets (`example*` objects). |
| `data-raw/` | Scripts that regenerate `data/` and some `inst/extdata/` files. |
| `inst/extdata/` | Small fixture files (GFF3, MAF, train) used by examples, and `simulatedRearrangements.rds` (200 simulated genome pairs stored as permutation vectors; see `?simulateRearrangements`). |
| `inst/CITATION` | Published-paper citation. |
| `man/` | roxygen-generated `.Rd` (do not edit by hand). |
| `vignettes/` | Long-form documentation (`*.Rmd`). |
| `CODE_REVIEW.md` | Reusability/code-review study. Excluded from the build. |

## Build / check / document

Generated artefacts (`man/`, `NAMESPACE`, `src/RcppExports.cpp`) are produced by
roxygen2/Rcpp — **edit the source, then regenerate**, never edit by hand.

```sh
R -e 'devtools::document()'      # regenerate man/ + NAMESPACE from roxygen
R -e 'devtools::load_all()'      # load for interactive work (compiles src/)
R CMD build .                    # build the source tarball
R CMD check GenomicBreaks_*.tar.gz   # full package check
R -e 'pkgdown::build_site()'     # rebuild the documentation site
```

Building the MAF reader requires a C++ toolchain (Rcpp). System prerequisites
for a clean environment are listed in `Singularity.def`.

GitHub Actions currently only builds/deploys the pkgdown site
(`.github/workflows/pkgdown.yaml`). It does not run `R CMD check`.

## Conventions

- **roxygen2 in markdown mode** (`Config/roxygen2`, `DESCRIPTION` `Roxygen:`).
  Document above each function; regenerate with `devtools::document()`.
- **`@family` / `@concept` tags matter**: the pkgdown reference index is grouped
  by them (`_pkgdown.yml`). Keep them on new exported functions.
- **S4 classes/methods** via `methods::setClass`/`setGeneric`/`setMethod`.
- **Semantic versioning** since 0.17.0; record every user-visible change in
  `NEWS.md` and bump the `Version:` in `DESCRIPTION`.
- Keep roxygen `@examples` runnable on the bundled example data — they are
  currently a primary form of verification.

## Structural-variant flags

### Conventions

- **Flags are set on the first range of a pattern** in _target_ order (for
  example the left flank of an inversion triplet).  The `filter*`, `remove*`,
  `show*` and `flip*` helpers find the other ranges by their offset from the
  flag.  This has confused users and may change later; keep it for now.
- **Look-ahead coding pattern** (see `flagDoubleInversions()`):
  - `qnext0 <- precede(gb$query) - seq_along(gb$query)` is the offset, in
    _target_ order, of the range that range `i` immediately precedes on the
    _query_; `qnextK <- c(tail(qnext0, -K), NA, ...)` is the same for range
    `i + K`, and `strandK` is built the same way from `strand(gb)`.
  - A pattern is an `&` chain of `qnextK == value` and `strandK == "+"` or
    `"-"` tests, piped to `sapply(isTRUE)` to turn `NA` into `FALSE`.
  - The expected `qnext` values of a pattern are
    `precede(IPos(r)) - seq_along(r)`, where `r` is the _query_ rank of each
    range in _target_ order.  The `NA` value (last range in _query_ order) is
    not tested.  Write this line as a comment above each pattern.
  - All ranges of a pattern must be on one _target_ sequence; the _query_
    adjacencies already imply one _query_ sequence.

### Plan (agreed 2026-10-23)

1. Fix the defects of the existing flagging functions with direct patches,
   one commit per issue.
2. Then replace the hand-written patterns by a table generated by
   enumeration, matched by one function (compare `allPossiblePairClasses` in
   `R/flagPairs.R`).  For `k` central ranges between two flanks that did not
   move relative to each other, enumerate all signed arrangements, keep those
   whose internal junctions are all breakpoints (coalesced input), and record
   for each its `qnextK`/`strandK` pattern, number of breakpoints `b`, minimal
   number of events `d` and breakpoint reuse rate `2d/b` (Pevzner and Tesler
   2003).  Each arrangement then gets exactly one label.
3. Parsimony: **a translocation counts as one event**, like an inversion.  An
   arrangement explained by fewer events with a translocation must not be
   labelled as several inversions.  Which moves count as one translocation
   (one range or a segment, inverted or not) is still to be decided.
4. First enumeration (3 central ranges, inversions only, breadth-first
   search; lower case means minus strand, flanks on the plus strand): of the
   48 arrangements, 25 have 4 breakpoints.  `b c A` and `C a b` need two
   overlapping inversions (`flagDoubleInversions()`); `a B c` is two separate
   inversions and `c B a` two nested ones (both found by `flagInversions()`);
   the other 21 need three inversions (reuse rate 1.5), and are flagged only
   as translocations (10) or not at all (11).  With 2 central ranges after
   coalescing: `a b` is a twin inversion, `b A` and `B a` need two inversions
   and `B A` three, and these three are flagged only as translocations.
   Triple inversions are to be revisited once translocations are included in
   the enumeration.
5. Deferred: a function returning, for each inversion, both breakpoint
   regions on both genomes, computed from the coordinates of the triplet.
   `leftInversionGaps()` only returns the left one on the _target_ genome.
6. The `show*`, `filter*` and `remove*` helpers must work well and
   consistently for all structural variants (see the known issues below).

### Known issues

**Flags are relative to the _target_ side, and fragile.**  A flag is a
metadata column of the _target_ ranges, and its meaning depends on the row
order of the object.  Any transformation of a `GBreaks` object (subsetting,
`sort()`, `c()`, `coalesce_contigs()`, `reverse()`, `swap()`) can invalidate
the flags or drop them, without warning.  `swap()` drops the _target_
metadata columns, flags included, and the `remove*` helpers drop them through
`coalesce_contigs()`.  Flags should be recomputed after a transformation.

**Translocations** (`flagTranslocations()`; tested against a brute-force
implementation of its definition on the stored simulations):

- The function finds transpositions (moves of one range, possibly inverted,
  on the same sequence or not) rather than translocations; do not rename it
  for now.
- A move of range `A` is seen in two places of the _target_ order, both
  flagged on the _target_ ranges, without `swap()`:
  - the **payload region** `U A V`: three consecutive _target_ ranges, with
    `U` and `V` adjacent and colinear on the _query_.  `tra` flag, on `U`;
    the payload is the next range.
  - the **site** `X Y`: two consecutive _target_ ranges between which `A` is
    found on the _query_ (`X A Y`).  `trs` flag, on `X`.  The pattern uses
    one indirect lookup: the payload is `i + qnext0[i]` (`i + qprev0[i]` on
    the minus strand), and the test is `qnext0[i] != 1 &
    qnext0[i] + qnext0[i + qnext0[i]] == 1`.
  When both are seen, they point to the same payload; either can be missing,
  for instance when `A` is alone on its sequence in one genome.  Events can
  be built by pairing the two flags through their payload.
- `flagAll()` only uses `tra`, on purpose for now: the `X` of a site can
  carry another flag, which would make its `stopifnot()` checks fail.
- Before 0.25.0, the site was flagged with `both = TRUE`, by a `tra` flag on
  the _query_ ranges computed through `swap()`; that flag was on `Y` for
  minus-strand sites, and on the wrong row when _target_ ranges were
  duplicated.
- `tol` limits the distances between the ranges that must be adjacent:
  `U A V` on the _target_ and `U V` on the _query_ for payload regions,
  `X Y` on the _target_ and `X A Y` on the _query_ for sites.
- The `show*`, `filter*` and `remove*` helpers use both flags; the moved
  ranges are computed by the internal `translocationPayloads()`.

**Strandless ranges** (strand `*`).  They come from (1) objects that are
unoriented as a whole, for instance made from reciprocal best hits, and (2)
_mapped regions_ (`bridgeRegions()`), which materialise the unaligned region
between two colinear alignments, possibly with a width of zero on one or both
genomes (in Plessy et al. 2024, some are new short introns).  Rule: process
when possible, ignore when needed, error when unavoidable.  Current state
(checked on unstranded examples and on examples with their bridge regions):

- Inversions (simple, double, twin) are never flagged when a range of the
  pattern is unstranded (done).
- Colinearity treats `*` as compatible with both strands, as
  `GenomicRanges` does, and adjacency in either direction counts for `*`
  pairs.  Works on unoriented objects and on mapped regions of non-zero
  width.
- `coalesce_contigs()` merges unoriented objects.  In mixed objects, it
  coalesces the oriented ranges and absorbs the mapped regions that they
  contain on both genomes (done; this avoids the zero-width problem below).
- Zero-width ranges: `precede()` reports a zero-width range as preceding
  itself, which breaks the adjacency tests of all flagging functions (a
  zero-width bridge is not colinear with the next alignment).
- `flagTranslocations()` treats `*` as `+`, so on an unoriented object the
  result depends on the arbitrary orientation of the _query_ sequences, and
  `flagAll()` can then fail.  With adjacency in either direction instead,
  extra ranges are flagged as moved (target `X Y U A V`, query `X A Y U V`:
  `Y` and `U` too, because moving `A` over `Y U` is also moving `Y U` over
  `A`); a stronger criterion is needed for unoriented data.
- `flagPairs()` returns `Bnd` (boundary) for every pair that involves `*`.
- `permutationVector()` counts `*` as `+`, so `inversionDistance()` returns
  a value for unoriented objects.
- To decide: whether mapped regions should be ignored by the flagging
  functions and `coalesce_contigs()` (computing on the stranded ranges only,
  which conflicts with flags located by row offsets), or make them stop with
  an explicit error; and how to flag translocations and pairs in unoriented
  objects.

**Helpers** (`show*`, `filter*`, `remove*`, `flip*`):

- Double inversions: there is no `showDoubleInversions()`.
- Twin inversions: there are no `show*`, `filter*` or `remove*` functions.
- `filterColinearRegions()` has no `detect` argument.
- The `detect` argument is documented as "reuse the flag or fail if it is
  absent", but the helpers return an empty object instead of failing; for
  the `remove*` functions, an empty object is especially misleading.
- Only `show*` and `filter*` have a `rename` argument.
- `removeDoubleInversions()` returns early for objects of fewer than 3
  ranges, instead of 5.

## Testing (please help here)

Verification today rests heavily on **vignettes + man-page examples**.

Unit tests are in `tests/testthat/` (**`testthat`**, edition 3); run them with
`R -e 'devtools::test()'`. Do **not** remove the existing vignette/@examples
coverage; add tests incrementally.

For structural variants, `simulateRearrangements()` and
`permutationVectorToGBreaks()` build test objects from known events or
signed permutations, and the simulations in
`inst/extdata/simulatedRearrangements.rds` can be compared with an
independent, brute-force implementation of a flag's definition.

Lift existing `@examples` into `expect_*` assertions. High-value candidates
include the colinearity/coalescing functions and the nucleotide-distance
functions (some distance functions are explicitly flagged as not-yet-proofread
in `NEWS.md`).
