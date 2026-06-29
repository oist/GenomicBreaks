# GenomicBreaks — Code Review / Reusability Study

Date: 2026-06-29
Package version reviewed: 0.22.0 (`DESCRIPTION:1-3`)
Reviewer: automated study (OpenCode), commissioned by the maintainer.

> **Status note.** There is **no official public release** of this package yet.
> Most engineering findings below are therefore **low urgency** and meant to be
> fixed opportunistically, not as release blockers. This document exists to
> *keep track* of them so they are not forgotten. Some entries may later be
> converted into GitHub issues.

---

## 1. Overview

`GenomicBreaks` is an R package that analyses **pairwise whole-genome
alignments** in which gene order has been scrambled by evolution. It is built
on the Bioconductor ecosystem (`GenomicRanges`, `Biostrings`, `BSgenome`,
`rtracklayer`, `S4Vectors`, `IRanges`, `GenomeInfoDb`, `BiocGenerics`,
`BiocParallel`) but is **not itself part of Bioconductor**. It is designed to
consume the output of nf-core's `pairgenomealign` pipeline, but can import
other pipelines' output as well (`README.md:15-17`).

Citation target: Plessy *et al.*, *Genome Research* 2024,
doi:10.1101/gr.278295.123 (`inst/CITATION`).

## 2. Architecture

### The `GBreaks` class

The core abstraction is the S4 class `GBreaks`, a thin wrapper around
`GenomicRanges::GRanges` (`R/GenomicBreaks.R:22`):

- The **main ranges** hold the coordinates on the **target** genome.
- A `query` **metadata column** (itself a `GRanges`) holds the matching
  coordinates on the **query** genome.
- Strand information is carried on the target ranges; query ranges are kept on
  the plus strand by convention (`R/load_genomic_breaks.R:16-17`).

A custom `initialize` method lets a `GBreaks` object be built directly from
`target` / `query` / `strand` `GRanges` (`R/GenomicBreaks.R:24-36`).

### Data flow

1. **Load** — `load_genomic_breaks()` dispatches on file extension to a GFF3
   reader (`rtracklayer::import.gff3`) or a MAF reader
   (`R/load_genomic_breaks.R:55-74`). The MAF reader is implemented in C++ via
   Rcpp (`src/readMaf.cpp`, `src/RcppExports.cpp`), with code partly derived
   from LAST (M. Frith, GPLv3+; see `NEWS.md:106-107`).
2. **Coalesce** — `coalesce_contigs()` / `flagColinearAlignments()` merge
   colinear alignment fragments (`R/coalesce_contigs.R`).
3. **Flag / filter** — structural variants: inversions, translocations, double
   and twin inversions (`flag*`, `filter*`, `show*`, `remove*` families;
   `NAMESPACE:30-78`).
4. **Quantify** — synteny / correlation / karyotype / tau indices and a large
   family of nucleotide distances (JC69, K80, K80+gap, F81, T92, HKY85, TN93,
   logDet, P-distance; `DESCRIPTION:50-119`).
5. **Visualise** — Oxford/dot plots and pair-of-chromosome plots via
   `genoPlotR` and `ggplot2` (`R/makeOxfordPlots.R`, `R/genoPlotR.R`).

### Documentation & distribution

- roxygen2 (markdown mode) + pkgdown site at <https://oist.github.io/GenomicBreaks>
  (`_pkgdown.yml:1`, `DESCRIPTION:40`).
- The pkgdown reference index is organised by roxygen `@family`/`@concept` tags
  (`_pkgdown.yml:5-72`). Those tags are functionally important, not cosmetic.
- Lazy-loaded example datasets in `data/` (`DESCRIPTION:17`).
- A Singularity recipe for reproducible installation (`Singularity.def`).

## 3. Findings

Severity legend: **High** (correctness/legal), **Medium** (reproducibility,
maintainability), **Low** (polish). Given the pre-release status, none of these
are urgent.

| # | Severity | Finding | Evidence | Suggested remediation |
|---|----------|---------|----------|-----------------------|
| 1 | Medium | **No automated test suite.** There is no `tests/` directory and no test framework configured. Verification relies entirely on vignettes and roxygen `@examples`. | no `tests/`; `vignettes/`, `man/` only | Add a `tests/testthat/` suite (see §5). |
| 2 | High | **Invalid `LICENSE` file.** Contains placeholder text instead of a real license/year/holder, although `DESCRIPTION` declares `BSD_3_clause + file LICENSE`. | `LICENSE:1-3` (`"CW and CP, or maybe just OIST?"`) | Decide the copyright holder and write a proper BSD-3-Clause `LICENSE` with year and holder. |
| 3 | Medium | **CI only builds the pkgdown site; no `R CMD check`.** Regressions in code/examples are not caught automatically. | `.github/workflows/` contains only `pkgdown.yaml` | Add an `R-CMD-check` GitHub Action (e.g. `r-lib/actions`), ideally on a Bioconductor-devel container. |
| 4 | Medium | **Unpinned dependencies.** No minimum versions on Bioconductor imports; installation pulls "latest", so results may drift across Bioc releases. | `DESCRIPTION:18-48`; `Singularity.def:24-25` | Record tested Bioc release; optionally add `>=` floors for packages with known API needs. |
| 5 | Medium | **Unproofread code shipped.** `TN93_distance()`, `HKY85_distance()`, `logDet_distance()` are documented as not yet proofread. | `NEWS.md:14-16` | Add validation tests against known references; mark experimental in docs until verified. |
| 6 | Medium | **No archival DOI; installation tracks a moving branch.** Users install HEAD of `main` via `install_github`, so a given analysis cannot be reproduced by version. | `README.md:29` | Tag releases (already semver since 0.17.0) and archive tagged versions (e.g. Zenodo) to mint a DOI. |
| 7 | Medium | **Non-reproducible container recipe.** `Bootstrap: docker / From: debian:trixie` plus `install_github(... main ...)` resolves differently over time. | `Singularity.def:1-2,25` | Pin the base image digest and install a tagged GenomicBreaks release; consider pinning R package versions. |
| 8 | Low | **roxygen stubs in the public API.** Several exported methods have `@return tbd` and `@param ... etc`. | `R/GenomicBreaks.R:113,150-156` | Fill in `@param`/`@return` for exported functions; an `R CMD check` (finding #3) would surface documentation gaps. |

## 4. Strengths to preserve

When addressing the above, take care **not** to regress these:

- Clean S4 design layered on a stable Bioconductor foundation.
- Rich, runnable roxygen `@examples` on lazy-loaded example data.
- A real `inst/CITATION` pointing at the published paper.
- Semantic versioning with a well-maintained `NEWS.md` (since 0.17.0).
- A container recipe and clear `README` install instructions.
- pkgdown reference grouping driven by `@family`/`@concept` tags.

## 5. Testing recommendation

The package currently leans on **vignettes + man-page examples** for
verification. This is valuable coverage and should be **kept**. In addition, the
recommendation is to **incrementally add a unit-test suite** using
**`testthat`**, which is the most frequently used test framework in Bioconductor
(the older `RUnit` convention is the legacy alternative; `tinytest` is a lighter-weight option).

Suggested approach:

1. Scaffold the standard layout:
   - `tests/testthat.R` — runner calling `testthat::test_check("GenomicBreaks")`.
   - `tests/testthat/test-*.R` — one file per topic.
   - Add `testthat` to `Suggests` in `DESCRIPTION`.
2. **Seed tests by lifting existing `@examples` into assertions.** Many examples
   already exercise representative behaviour and can become `expect_*` checks
   with little effort. Good first candidates:
   - `isSorted()` / `range()` / `subsetByOverlaps()` on the example data
     (`R/GenomicBreaks.R:83-84,119-120,162`).
   - `coalesce_contigs()` / `flagColinearAlignments()` colinearity behaviour,
     including the minus-strand and strandless cases already demonstrated
     (`R/coalesce_contigs.R:32-58`).
   - `load_genomic_breaks()` round-trips from the bundled `inst/extdata`
     fixtures (GFF3 and MAF), guarding against regressions like the historical
     1-nt MAF coordinate shifts (`NEWS.md:91-93,131-133`).
   - Distance functions against hand-computed or literature reference values —
     this would also address finding #5 for the unproofread distances.
3. Wire the suite into CI alongside the `R-CMD-check` action (finding #3).

This complements, rather than replaces, the example/vignette coverage: examples
keep the docs trustworthy, while `tests/testthat/` pins down edge cases and
numerical results.

## 6. Pointers

- Orientation for future contributors / AI sessions: `AGENTS.md`.
- This document is excluded from the package tarball via `.Rbuildignore`
  (`^CODE_REVIEW\.md$`).
