# AGENTS.md — Orientation for GenomicBreaks

This file orients contributors and AI coding sessions working on
`GenomicBreaks`. For a detailed review of known engineering gaps, see
[`CODE_REVIEW.md`](CODE_REVIEW.md).

## Status / urgency

**There is no official public release of this package yet.** Consequently, the
engineering issues recorded in `docs/code_review.md` (license placeholder,
missing CI, unpinned dependencies, etc.) are **low urgency**: fix them
opportunistically, not as blockers. Do not let them stall feature work.

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
| `data-raw/` | Scripts that regenerate `data/`. |
| `inst/extdata/` | Small fixture files (GFF3, MAF, train) used by examples. |
| `inst/CITATION` | Published-paper citation. |
| `man/` | roxygen-generated `.Rd` (do not edit by hand). |
| `vignettes/` | Long-form documentation (`*.Rmd`). |
| `docs/code_review.md` | Reusability/code-review study. Excluded from the build. |

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

## Testing (please help here)

Verification today rests heavily on **vignettes + man-page examples**. The
maintainer is keen to grow real unit tests. When touching code, **proactively
suggest adding or migrating coverage into a `tests/testthat/` suite** using
**`testthat`** — the most common test framework in Bioconductor (`RUnit` is the
legacy alternative). Do **not** remove the example/vignette coverage; the test
suite should complement it.

A practical starting point is to lift existing `@examples` (e.g.
`coalesce_contigs()`, `isSorted()`, `load_genomic_breaks()` round-trips, the
distance functions) into `expect_*` assertions. See `docs/code_review.md` §5 for
concrete seed-test ideas. Note that the distance functions `TN93_distance()`,
`HKY85_distance()`, and `logDet_distance()` are flagged as not-yet-proofread
(`NEWS.md`), so reference-value tests there are especially valuable.
