# Changelog

## GenomicBreaks 0.26.0

- `GenomicBreaks` now imports `Seqinfo` instead of `GenomeInfoDb`, and
  therefore requires Bioconductor 3.22 or later. The example objects
  need the `Seqinfo` package, where the `Seqinfo` class moved in
  Bioconductor 3.22, but installation on earlier releases did not stop
  for lack of it
  ([\#70](https://github.com/oist/GenomicBreaks/issues/70)). The
  functions that `GenomicBreaks` took from `GenomeInfoDb` are all
  defined in `Seqinfo`, which `GenomeInfoDb` re-exports only
  temporarily.
- New *monkey indices*:
  [`wukong_index()`](https://oist.github.io/GenomicBreaks/reference/monkey_indices.md),
  [`abu_index()`](https://oist.github.io/GenomicBreaks/reference/monkey_indices.md),
  [`albert_index()`](https://oist.github.io/GenomicBreaks/reference/monkey_indices.md),
  [`marcel_index()`](https://oist.github.io/GenomicBreaks/reference/monkey_indices.md),
  [`imo_index()`](https://oist.github.io/GenomicBreaks/reference/monkey_indices.md),
  [`rafiki_index()`](https://oist.github.io/GenomicBreaks/reference/monkey_indices.md),
  [`hanuman_index()`](https://oist.github.io/GenomicBreaks/reference/monkey_indices.md),
  [`jack_index()`](https://oist.github.io/GenomicBreaks/reference/monkey_indices.md),
  [`able_index()`](https://oist.github.io/GenomicBreaks/reference/monkey_indices.md),
  [`mizaru_index()`](https://oist.github.io/GenomicBreaks/reference/monkey_indices.md),
  [`george_index()`](https://oist.github.io/GenomicBreaks/reference/monkey_indices.md)
  and
  [`baker_index()`](https://oist.github.io/GenomicBreaks/reference/monkey_indices.md).
  They are twelve summary statistics that were not designed to represent
  a biological fact: they were retained from a screen of 5,420
  mechanically assembled functions because their values are
  reproducible, structured by phylogeny, different from those of
  shuffled alignments, and not strongly correlated with the existing
  statistics nor with each other. See
  [`?monkey_indices`](https://oist.github.io/GenomicBreaks/reference/monkey_indices.md).
- New vignette *Monkey indices*
  ([`vignette("MonkeyIndices")`](https://oist.github.io/GenomicBreaks/articles/MonkeyIndices.md)),
  which defines the monkey indices and reports how they were selected
  and how they behave on the alignments of the screen. The code of the
  screen is in `inst/scripts/monkey_screen/`.
- [`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md)
  cuts the alignments of all windows in one call instead of one call per
  window (13 s to 0.2 s on the yeast example alignment).
- **Breaking change:**
  [`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md)
  now returns a `GRangesList`, as it already did with `cut = FALSE`,
  instead of a list of `GBreaks` objects.
- Fix
  [`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md),
  which did not cut the *query* ranges of objects loaded by
  [`load_genomic_breaks()`](https://oist.github.io/GenomicBreaks/reference/load_genomic_breaks.md):
  it read the orientation of the alignments on the *query* ranges
  instead of the *target* ranges.
- [`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md)
  accepts objects without a `score` column, and its example runs on the
  lazy-loaded `exampleInversion` object.

## GenomicBreaks 0.25.0

- The handling of ranges on the `*` strand (unoriented objects, and
  mapped regions such as the output of
  [`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md))
  is experimental: it may change in later versions without being
  considered a breaking change. Inversions are never flagged on such
  ranges, and functions that need orientation, such as
  [`permutationVector()`](https://oist.github.io/GenomicBreaks/reference/permutationVector.md),
  stop with an error.
- New
  [`simulateRearrangements()`](https://oist.github.io/GenomicBreaks/reference/simulateRearrangements.md)
  function, which simulates a pair of genomes rearranged by random
  inversions and transpositions and records the events, and new
  [`permutationVectorToGBreaks()`](https://oist.github.io/GenomicBreaks/reference/permutationVectorToGBreaks.md)
  function, the inverse of
  [`permutationVector()`](https://oist.github.io/GenomicBreaks/reference/permutationVector.md).
  A set of 200 simulations is provided in
  `inst/extdata/simulatedRearrangements.rds` for tests and benchmarks.
- **Breaking change:**
  [`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md)
  now flags the two places where a moved range is seen in the order of
  the *target* genome: its payload region (`tra` flag, as before) and
  its site, the pair of adjacent ranges between which it is found on the
  *query* genome (new `trs` flag, on the first range of the pair). The
  `both` argument and the `tra` flag on the *query* ranges are removed:
  they computed the site by swapping the genomes and sorting again,
  which put the flag on the second range of the site when it was on the
  minus strand, and on the wrong row when *target* ranges were
  duplicated.
  [`flagAll()`](https://oist.github.io/GenomicBreaks/reference/flagAll.md)
  still uses only the `tra` flag.
- The `tol` argument of
  [`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md)
  is now implemented, with the same meaning as in
  [`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md).
- [`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md)
  now coalesces over *mapped regions*: in an object that mixes oriented
  ranges and ranges on the `*` strand, such as the output of
  [`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md)
  added to the alignments, the oriented ranges are coalesced and absorb
  the unstranded ranges that they contain on both genomes
  (experimental). Before, such objects were not coalesced at all, and
  mapped regions of width zero on one genome caused an error.
- [`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md)
  no longer flags triplets that contain a range on the `*` strand, such
  as a colinear `+ * +` triplet made by adding the output of
  [`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md)
  to an object: inversions can not be determined without orientation.
  [`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md)
  and
  [`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md)
  already required known strands; this is now documented.
- [`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md)
  now tests unstranded ranges in both directions, instead of as
  plus-strand ranges, so that its result on unoriented objects no longer
  depends on the orientation of the *query* sequences. As orientation is
  unknown, more ranges can be flagged as moved (see the help page).
  [`flagAll()`](https://oist.github.io/GenomicBreaks/reference/flagAll.md)
  gives precedence to colinearity on such ranges, instead of failing
  (experimental).
- [`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md)
  no longer flags all pairs that involve a range on the `*` strand as
  boundaries (`Bnd`): they are now `Col` when adjacent in an order
  compatible with colinearity, `Scr` when not adjacent, and `NA` when
  adjacent but undetermined (experimental). The factor levels are
  unchanged.
- Ranges of width zero on one or both genomes, such as the mapped
  regions that
  [`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md)
  returns between adjacent alignments, broke the adjacency tests of
  [`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md),
  [`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md)
  and the functions that flag structural variants:
  [`GenomicRanges::precede()`](https://rdrr.io/pkg/IRanges/man/nearest-methods.html)
  and [`follow()`](https://rdrr.io/pkg/IRanges/man/nearest-methods.html)
  report a zero-width range as adjacent to itself. The adjacency tests
  now place such ranges between the bases that surround them. Results on
  objects without zero-width ranges are unchanged.
- [`permutationVector()`](https://oist.github.io/GenomicBreaks/reference/permutationVector.md),
  and therefore
  [`inversionDistance()`](https://oist.github.io/GenomicBreaks/reference/inversionDistance.md)
  and
  [`breakpointGraphProperties()`](https://oist.github.io/GenomicBreaks/reference/breakpointGraphProperties.md),
  now stop with an error when a range is on the `*` strand, instead of
  counting it as a plus-strand range.
- [`showTranslocations()`](https://oist.github.io/GenomicBreaks/reference/showTranslocations.md),
  [`filterTranslocations()`](https://oist.github.io/GenomicBreaks/reference/filterTranslocations.md)
  and
  [`removeTranslocations()`](https://oist.github.io/GenomicBreaks/reference/removeTranslocations.md)
  now use both the `tra` and `trs` flags, so that moves seen only at
  their site are no longer ignored. A moved range seen at both places is
  counted once.
- **Breaking change:**
  [`dist2next()`](https://oist.github.io/GenomicBreaks/reference/dist2next.md)
  now returns the number of bases between two ranges, as
  [`GenomicRanges::distance()`](https://rdrr.io/pkg/IRanges/man/nearest-methods.html)
  does (`0` for adjacent or overlapping ranges), instead of that number
  plus one. The `tol` argument of
  [`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
  [`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md)
  and
  [`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md)
  now bridges gaps of up to `tol` bases, as documented, instead of gaps
  shorter than `tol`. With a finite `tol`, results change where a gap is
  exactly `tol` bases wide; for example `coalesce_contigs(gb, tol = 0)`
  now merges adjacent colinear ranges. Results with the default
  `tol = Inf` are unchanged.
- [`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md)
  now requires the ranges before and after the triplet to be colinear
  with it: on the plus strand for the `bcA` and `Cab` patterns, and on
  the minus strand for `BAc` and `aCB`. Previously, it also flagged
  triplets that require three inversions, and triplets whose neighbours
  map elsewhere in the *query* genome. The flag is now set on the left
  flank, as in
  [`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md)
  and
  [`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md).
  As a consequence,
  [`filterDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/filterDoubleInversions.md)
  now returns the central triplet (it returned its last two ranges and
  the right flank), and
  [`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md)
  was adjusted to remove the same three ranges.
- [`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md)
  and
  [`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md)
  now require all the ranges of a pattern to be on the same *target*
  sequence, like
  [`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md)
  does.
- [`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md)
  and
  [`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md)
  ignored the first range of the object when testing whether two ranges
  are adjacent on the *query* genome. For example,
  `flagPairs(exampleInversionGarg2019)` now returns `Scr` instead of
  `Inv` for ranges 4 and 5, which are separated by range 1 on the
  *query*.
- `filterTranslocations(remove = TRUE)` returned an empty object when no
  translocation was flagged; it now returns all the ranges.
- [`leftInversionGaps()`](https://oist.github.io/GenomicBreaks/reference/leftInversionGaps.md)
  now returns one gap per inversion, as a zero-width range when the left
  flank and the central range are adjacent or overlap (these inversions
  were dropped), with a new `flagIndex` metadata column giving the
  position of the left flank. A new `detect` argument allows to reuse
  existing `inv` flags, for instance computed with a `tol` value.

## GenomicBreaks 0.24.0

- Support correction for unaligned bases in
  [`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md).
- Refactor breakpoint graph code.
- Add new example `GBreaks` objects from the litterature.

## GenomicBreaks 0.23.0

- Support for the *blasttab+* format in
  [`load_genomic_breaks()`](https://oist.github.io/GenomicBreaks/reference/load_genomic_breaks.md).

## GenomicBreaks 0.22.0

- New
  [`K80_gap_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_gap_distance.md)
  function for Nishimaki and Sato’s K2P + Gap distance.

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
- New
  [`subsetByOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)
  range that combines regions overlapping by their *target* or by their
  *query* ranges.

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
