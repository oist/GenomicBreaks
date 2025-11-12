# Changelog

## GenomicBreaks 0.21.0

- New
  [`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md)
  function with example data.
- Support all-zero substitution matrices in the distance functions
  introduced earlier.

## GenomicBreaks 0.20.0

- Correct GC content calculation in
  [`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md)
  and add two alternatives.
- Add three more distances suggested by Copilot. The code is not
  proofread yet; take the results with skepticism.
  [`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
  [`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
  [`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md).

## GenomicBreaks 0.19.0

- Add a
  [`GCequilibrium()`](https://oist.github.io/GenomicBreaks/reference/GCequilibrium.md),
  [`GCpressure()`](https://oist.github.io/GenomicBreaks/reference/GCpressure.md)
  and
  [`GCproportion()`](https://oist.github.io/GenomicBreaks/reference/GCproportion.md)
  function.
- Add a
  [`P_error()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md)
  function to compute the binomial standard error on
  [`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md)
- Add lazy-loaded example data, `exampleProbabilityMatrix` and
  `exampleSubstitutionMatrix` for easy testing of the functions above.

## GenomicBreaks 0.18.0

- New
  [`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)
  function. Closes
  [\#16](https://github.com/oist/GenomicBreaks/issues/16)
- Reform the nucleotide distance functions to also take count matrices
  as input.
- Implement the 4 variants of percent nucleotide difference in
  [`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md)
  as defined by [May (2004)](https://doi.org/10.1016/j.str.2004.04.001).
  Closes [\#26](https://github.com/oist/GenomicBreaks/issues/26)
- Add a
  [`gapProportion()`](https://oist.github.io/GenomicBreaks/reference/gapProportion.md)
  function.

## GenomicBreaks 0.17.0

- Switch to semantic versioning.
- New
  [`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md)
  (ungapped percent difference) distance function.

## GenomicBreaks 0.16.3

- New
  [`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md)
  (Felsenstein 1981) distance function.
- Documentation and namespace fixes.

## GenomicBreaks 0.16.2

- New
  [`matchPairs()`](https://oist.github.io/GenomicBreaks/reference/matchPairs.md)
  function to remove non-syntenic regions.

## GenomicBreaks 0.16.1

- Fix a bug where
  [`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
  [`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
  and
  [`removeTranslocations()`](https://oist.github.io/GenomicBreaks/reference/removeTranslocations.md)
  would discard all the ranges if no inversion, double inversion or
  translocation (respectively) were found.

## GenomicBreaks 0.16.0

- New functions to consolidate alignments by removing structural
  variants and re-coalescing the objects:
  [`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
  [`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
  and
  [`removeTranslocations()`](https://oist.github.io/GenomicBreaks/reference/removeTranslocations.md).

## GenomicBreaks 0.15.0

- New `detect` option in the show and filter functions for inversions
  and translocations, to re-compute flags automatically and simplify
  interactive use. Set to `FALSE` and flag manually when performance is
  an issue.

## GenomicBreaks 0.14.10

- New
  [`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md)
  function to compute similarity indices on windowed intervals.
- Experimental `tile` argument in
  [`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md).
  It will be removed once
  [`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md)
  can be applied to `strand_randomisation_index`.

## GenomicBreaks 0.14.9

- New distance functions written by Zikun Yang.

## GenomicBreaks 0.14.8

- Update for Bioc 3.21 and R 4.5.
- New experimental inversion analysis functions written by Bruna
  Fistarol.
- Namespace and documentation fixes.
- Guard against empty `GBreaks` objects in
  [`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md).

### Dropped features

- Removed `bp_coverage`, not used anymore.
- Removed support for Axt files as CNEr is broken for the moment.

## GenomicBreaks 0.14.7

- Correct a bug in
  [`load_genomic_breaks()`](https://oist.github.io/GenomicBreaks/reference/load_genomic_breaks.md)
  where sequence lengths would be too long of 1 nucleotide when loaded
  from MAF files.

## GenomicBreaks 0.14.6

- Experimental support for loading sequence information from FASTA
  files.

## GenomicBreaks 0.14.5

- Support loading MAF files where alignment sign is on the query.
- Correct example data to make it really one-to-one.

## GenomicBreaks 0.14.4

- Remove dependency on BOOST thanks to C++ code taken from LAST (M.
  Frith / GPLv3+).

## GenomicBreaks 0.14.3

### New functionalities

- New
  [`scaffoldByFlipAndMerge()`](https://oist.github.io/GenomicBreaks/reference/scaffoldByFlipAndMerge.md)
  and
  [`splitSeqLevel()`](https://oist.github.io/GenomicBreaks/reference/splitSeqLevel.md)
  functions for trivial scaffolding by hand.

## GenomicBreaks 0.14.2

### New functionalities

- New
  [`GBreaksToMatrix()`](https://oist.github.io/GenomicBreaks/reference/GBreaksToMatrix.md)
  function to apply image processing approaches on scrambling patterns.
- New
  [`keepLongestPair()`](https://oist.github.io/GenomicBreaks/reference/keepLongestPair.md)
  function to focus on the longest sequence feature in each genome.
- New
  [`wholeGenomeClassification()`](https://oist.github.io/GenomicBreaks/reference/wholeGenomeClassification.md)
  function to classify *isolated alignments*, *collinear alignments*,
  *bridge regions* and *breakpoint regions*.

## GenomicBreaks 0.14.1

### Bug fixes

- Correct a 1-nt shift in minus-strand coordinates of the *query* genome
  of *GBreaks* objects loaded from MAF files. This shift could cause
  overlapping ranges that crashed the coalescing algorithm.
- [`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md)
  now returns empty *GBreaks* objects if no bridge region was found.
  Closes issue [\#20](https://github.com/oist/GenomicBreaks/issues/20).
- Allow
  [`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md)
  to work on genomes larger than the maximum 32-bit integer. Closes
  issue [\#24](https://github.com/oist/GenomicBreaks/issues/24).

### New functionalities

- New
  [`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md)
  function to remove detected inversions before searching for nested
  ones.
- New
  [`filterDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/filterDoubleInversions.md)
  function for extracting double inversions.

## GenomicBreaks 0.14.0

- The
  [`load_genomic_breaks()`](https://oist.github.io/GenomicBreaks/reference/load_genomic_breaks.md)
  function gains support for files in the *Multiple Alignment Format*
  (MAF).

## GenomicBreaks 0.13.1

### New functionalities

- The
  [`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md)
  function gains a `col` argument for choosing to color by strand or by
  score.

### Bug fixes

- In
  [`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md),
  the sequence names are now centered on their region.

## GenomicBreaks 0.13.0

### Important bug fixes

- The
  [`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md)
  and
  [`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md)
  functions were corrected so that strandless ranges are called
  *colinear* when they satisfy the conditions on any strand. Before,
  they would be treated like plus-strand ranges, causing failure to
  detect colinearity on the minus strand.

### New functionalities

- New
  [`getSeq()`](https://oist.github.io/GenomicBreaks/reference/getSeq.md)
  method for `GRanges` and `GBreaks` objects that have proper `seqinfo`
  for an installed `BSgenome` package.
- New
  [`pairwiseAlignment()`](https://oist.github.io/GenomicBreaks/reference/pairwiseAlignment.md)
  method for `GRanges` and `GBreaks` objects, that will call
  [`getSeq()`](https://oist.github.io/GenomicBreaks/reference/getSeq.md)
  automatically.
- New “mid” option to the `direction` argument of
  [`get_bps()`](https://oist.github.io/GenomicBreaks/reference/get_bps.md),
  to output breakpoints at mid-distance between ranges.
- New
  [`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md)
  function that maps unaligned regions of the *target* to the *query*
  genome when they are flanked by colinear regions.
- [`filterTranslocations()`](https://oist.github.io/GenomicBreaks/reference/filterTranslocations.md)
  gains a `remove` flag and `flagTranslocations` gains a `both` flag.
- `swap` gains a `sort` flag, so that in interactive sessions,
  `gb |> swap() |> sort(i=T)` can be replaced by the shorter
  `gb |> swap(s = T)`.

### Backwards-incompatible changes

- The function `bp_heatmap` now expects to genome sequences to be
  reachable via
  [`BSgenome::getBSgenome`](https://rdrr.io/pkg/BSgenome/man/available.genomes.html).
- The function `dist2next` now operates on `GRanges` and `GBreaks`, and
  its main argument is renamed to `x`.

## GenomicBreaks 0.12.2

### New functionalities

- The function `forceSeqLevels()` now operates on both the *target* and
  the *query* ranges of `GBreaks` objects.
- New example data set: *Neurospora crassa* chromosome III aligned to
  its homologue in *Podospora comata*: chromosome 7.
- New
  [`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md).
- [`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md)
  includes strand information in the returned object.
- New `type = 'none'` option for
  [`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md),
  for better customisation.

## GenomicBreaks 0.12.1

- The
  [`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)
  function now plots the name of the selected sequence levels instead of
  `target` and `query`. The labels can be overriden with the
  `dna_seg_labels` parameters. Closes issue
  [\#6](https://github.com/oist/GenomicBreaks/issues/6).

## GenomicBreaks 0.12.0

### Backwards-incompatible changes

- Removed `selChroms` parameter in
  [`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md).
- New defaults for
  [`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md):
  fixed-ratio line plot with round ends and size of 1. Axis labels are
  now “target” and “query” by default. In any case, Seqname is displayed
  instead of genome name if there is only one sequence in the ranges
  object.

### New functionalities

- New function
  [`orderQuerySeqLevels()`](https://oist.github.io/GenomicBreaks/reference/orderQuerySeqLevels.md)
  to better diagonalise dot plots.
- New `diag` parameter in
  [`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md)
  to diagonalise the plots by re-ordering *query* sequence levels.
- New `size` parameter in
  [`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md)
- New
  [`guessSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/guessSeqLengths.md)
  and
  [`forceSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/forceSeqLengths.md)
  functions to use as last ressort when sequence lengths are not known.

## GenomicBreaks 0.11.0

### Backwards-incompatible changes

- Remove project-specific `load_one_to_ones` function.

## GenomicBreaks 0.10.0

### Backwards-incompatible changes

- Removed `flagInsersions`, redundant with
  [`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md)
  and buggy on minus strands.
- [`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)
  now plots coordinates and a scale bar by default.

### New functionalities

- New
  [`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md)
  function to exhaustively flag windows of 2 successive `GBreaks`
  entries.
- New
  [`gb2xlim()`](https://oist.github.io/GenomicBreaks/reference/gb2xlim.md)
  function to narrow focus on plots made with
  [`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md).
- New
  [`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md)
  function to create macrosynteny plots. Thanks to Aleksandra Bliznina.
- In
  [`mergeSeqLevels()`](https://oist.github.io/GenomicBreaks/reference/mergeSeqLevels.md),
  use maximal end value for each `seqlevel` if no if no `seqlengths`
  were provided.
- A new set of example 5-uplets, to illustrate double inversions.
- Flag double inversions with
  [`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md)
  and
  [`flagAll()`](https://oist.github.io/GenomicBreaks/reference/flagAll.md),
  test them with `exampleDoubleInversion`.

### Bug fixes

- Add strand check on
  [`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md)
  so that it does not return a false positive on
  `exampleDoubleInversion`.

## GenomicBreaks 0.9.1

### New functionalities

- [`flagAll()`](https://oist.github.io/GenomicBreaks/reference/flagAll.md)
  now ensures that flags do not collide, and reports translocations.

## GenomicBreaks 0.9.0

### New functionalities

- New functions
  [`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
  [`showTranslocations()`](https://oist.github.io/GenomicBreaks/reference/showTranslocations.md)
  a
  [`filterTranslocations()`](https://oist.github.io/GenomicBreaks/reference/filterTranslocations.md).
- New
  [`filterColinearRegions()`](https://oist.github.io/GenomicBreaks/reference/filterColinearRegions.md)
  function.
- New
  [`cleanGaps()`](https://oist.github.io/GenomicBreaks/reference/cleanGaps.md)
  function that is guaranteed to return only the intervals between the
  input ranges.
- The
  [`GBreaks()`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  constructor now accepts three new parameters, `target`, `query` and
  `strand`, to construct `GBreaks` objects directly from `GRanges`
  objects.
- New [`range()`](https://rdrr.io/r/base/range.html) methods for
  `GBreaks` objects, that runs `GRanges::range` on each pair of
  target/query *seqnames*, ignoring strand information.
- New `subsetByOverlaps()` range that combines regions overlapping by
  their *target* or by their *query* ranges.

### Bug fixes

- [`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md)
  now handles properly strandless ranges. (Issue
  [\#1](https://github.com/oist/GenomicBreaks/issues/1))
- [`swap()`](https://oist.github.io/GenomicBreaks/reference/swap.md) now
  returns proper
  [`GBreaks()`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  objects.

## GenomicBreaks 0.8.0

### Backwards-incompatible changes

- Remove `flagInterruptions`, which was nonsense.

### New functionalities

- New
  [`isSorted()`](https://oist.github.io/GenomicBreaks/reference/isSorted.md)
  method to properly test with `ignore.strand` set to `TRUE`.
- Default
  [`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)
  to the first *sequence level* of the `GBreaks` object.
- [`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)
  now renders regions as blocks and not arrows as the strand information
  is encoded in the color and shape of the comparison bands. It also
  plots in light blue the blocks that are not matched. Closes issues
  [\#3](https://github.com/oist/GenomicBreaks/issues/3) and
  [\#5](https://github.com/oist/GenomicBreaks/issues/5).
- Added example data, see ?StructuralVariants and \`vignette(“Structural
  Variants”, package = “GenomicBreaks”) for details.

### Bug fixes

- In
  [`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
  check if input was sorted. Closes Issue
  [\#2](https://github.com/oist/GenomicBreaks/issues/2).
- [`flagLongShort()`](https://oist.github.io/GenomicBreaks/reference/flagLongShort.md)
  now handles the case where some annotations do not match.

## GenomicBreaks 0.7.1

- Fix namespace and variable bugs that were preventing
  [`reverse()`](https://oist.github.io/GenomicBreaks/reference/reverse.md)
  to run.

## GenomicBreaks 0.7.0

- New
  [`mergeSeqLevels()`](https://oist.github.io/GenomicBreaks/reference/mergeSeqLevels.md)
  function to concatenate sequence features, for instance before
  plotting.
- New
  [`reverse()`](https://oist.github.io/GenomicBreaks/reference/reverse.md)
  method for `GBreaks` and `GRanges` objects.
- Really return a `GBreaks` object in
  [`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md).
- Pass extra arguments from
  [`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)
  to
  [`genoPlotR::plot_gene_map`](https://rdrr.io/pkg/genoPlotR/man/plot_gene_map.html).

## GenomicBreaks 0.6.0

- Add pairwise chromosome plot functions using the
  [`genoPlotR`](http://genoplotr.r-forge.r-project.org) package, see in
  particular
  [`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md).

## GenomicBreaks 0.5.3

- *pkgdown* update.

## GenomicBreaks 0.5.2

- Greatly speed up
  [`flagLongShort()`](https://oist.github.io/GenomicBreaks/reference/flagLongShort.md).

## GenomicBreaks 0.5.1

- Correct a variable name bug in
  [`flagLongShort()`](https://oist.github.io/GenomicBreaks/reference/flagLongShort.md).
- Export the
  [`GBreaks()`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  constructor.

## GenomicBreaks 0.5.0

- Add a
  [`flagLongShort()`](https://oist.github.io/GenomicBreaks/reference/flagLongShort.md)
  function to transfer annotations about chromosome arms.

## GenomicBreaks 0.4.0

- Use alignments of *Saccharomyces paradoxus* to *S. cerevisiae* as
  example data.
- Addition of a
  [`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md)
  function to calculate Gene Order Conservation.

## GenomicBreaks 0.3.0

- Added the
  [`leftInversionGaps()`](https://oist.github.io/GenomicBreaks/reference/leftInversionGaps.md)
  function.

## GenomicBreaks 0.2.0

- Document the package with [`pkgdown`](https://pkgdown.r-lib.org/).
