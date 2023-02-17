# GenomicBreaks Unreleased

## Bug fixes

 * `bridgeRegions()` now returns empty _GBreaks_ objects if no bridge region
   was found.  Closes issue #20.

## New functionalities

 * New `flipInversions()` function to remove detected inversions before
   searching for nested ones.

# GenomicBreaks 0.14.0

 * The `load_genomic_breaks()` function gains support for files in the
   _Multiple Alignment Format_ (MAF).

# GenomicBreaks 0.13.1

## New functionalities

 * The `makeOxfordPlots()` function gains a `col` argument for choosing to color
   by strand or by score.
   
## Bug fixes

 * In `makeOxfordPlots()`, the sequence names are now centered on their region.

# GenomicBreaks 0.13.0

## Important bug fixes

 * The `flagColinearAlignments()` and `coalesce_contigs()` functions were
   corrected so that strandless ranges are called _colinear_ when they satisfy
   the conditions on any strand.  Before, they would be treated like plus-strand
   ranges, causing failure to detect colinearity on the minus strand.

## New functionalities

 * New `getSeq()` method for `GRanges` and `GBreaks` objects that have proper
   `seqinfo` for an installed `BSgenome` package.
 * New `pairwiseAlignment()` method for `GRanges` and `GBreaks` objects, that
   will call `getSeq()` automatically.
 * New "mid" option to the `direction` argument of `get_bps()`, to output
   breakpoints at mid-distance between ranges.
 * New `bridgeRegions()` function that maps unaligned regions of the _target_ to
   the _query_ genome when they are flanked by colinear regions.
 * `filterTranslocations()` gains a `remove` flag and `flagTranslocations`
   gains a `both` flag.
 * `swap` gains a `sort` flag, so that in interactive sessions,
   `gb |> swap() |> sort(i=T)` can be replaced by the shorter `gb |> swap(s = T)`.

## Backwards-incompatible changes

 * The function `bp_heatmap` now expects to genome sequences to be reachable
   via `BSgenome::getBSgenome`.
 * The function `dist2next` now operates on `GRanges` and `GBreaks`, and its
   main argument is renamed to `x`.

# GenomicBreaks 0.12.2

## New functionalities

 * The function `forceSeqLevels()` now operates on both the _target_ and
   the _query_ ranges of `GBreaks` objects.
 * New example data set:  _Neurospora crassa_ chromosome III aligned to its
   homologue in _Podospora comata_: chromosome 7.
 * New `strand_randomisation_index()`.
 * `makeOxfordPlots()` includes strand information in the returned object.
 * New `type = 'none'` option for `makeOxfordPlots()`, for better customisation.

# GenomicBreaks 0.12.1

 * The `plotApairOfChrs()` function now plots the name of the selected sequence
   levels instead of `target` and `query`.  The labels can be overriden with
   the `dna_seg_labels` parameters.  Closes issue #6.

# GenomicBreaks 0.12.0

## Backwards-incompatible changes

 * Removed `selChroms` parameter in `makeOxfordPlots()`.
 * New defaults for `makeOxfordPlots()`: fixed-ratio line plot with round ends and
   size of 1.  Axis labels are now "target" and "query" by default.  In any
   case, Seqname is displayed instead of genome name if there is only one
   sequence in the ranges object.

## New functionalities

 * New function `orderQuerySeqLevels()` to better diagonalise dot plots.
 * New `diag` parameter in `makeOxfordPlots()` to diagonalise the plots by
   re-ordering _query_ sequence levels.
 * New `size` parameter in `makeOxfordPlots()`
 * New `guessSeqLengths()` and `forceSeqLengths()` functions to use as last
   ressort when sequence lengths are not known.

# GenomicBreaks 0.11.0

## Backwards-incompatible changes

* Remove project-specific `load_one_to_ones` function.

# GenomicBreaks 0.10.0

## Backwards-incompatible changes

* Removed `flagInsersions`, redundant with `flagTranslocations()` and buggy
  on minus strands.
* `plotApairOfChrs()` now plots coordinates and a scale bar by default.

## New functionalities

* New `flagPairs()` function to exhaustively flag windows of 2 successive
  `GBreaks` entries.
* New `gb2xlim()` function to narrow focus on plots made with `plotApairOfChrs()`.
* New `makeOxfordPlots()` function to create macrosynteny plots.  Thanks to
  Aleksandra Bliznina.
* In `mergeSeqLevels()`, use maximal end value for each `seqlevel` if no
  if no `seqlengths` were provided.
* A new set of example 5-uplets, to illustrate double inversions.
* Flag double inversions with `flagDoubleInversions()` and `flagAll()`, test
  them with `exampleDoubleInversion`.

## Bug fixes

* Add strand check on `flagTranslocations()` so that it does not return a
  false positive on `exampleDoubleInversion`.

# GenomicBreaks 0.9.1

## New functionalities

* `flagAll()` now ensures that flags do not collide, and reports translocations.

# GenomicBreaks 0.9.0

## New functionalities

* New functions `flagTranslocations()`, `showTranslocations()` a
  `filterTranslocations()`.
* New `filterColinearRegions()` function.
* New `cleanGaps()` function that is guaranteed to return only the intervals
  between the input ranges.
* The `GBreaks()` constructor now accepts three new parameters, `target`,
  `query` and `strand`, to construct `GBreaks` objects directly from `GRanges`
  objects.
* New `range()` methods for `GBreaks` objects, that runs `GRanges::range` on
  each pair of target/query _seqnames_, ignoring strand information.
* New `subsetByOverlaps()` range that combines regions overlapping by their
  _target_ or by their  _query_ ranges.

## Bug fixes

* `flagColinearAlignments()` now handles properly strandless ranges. (Issue #1)
* `swap()` now returns proper `GBreaks()` objects.

# GenomicBreaks 0.8.0

## Backwards-incompatible changes

* Remove `flagInterruptions`, which was nonsense.

## New functionalities

* New `isSorted()` method to properly test with `ignore.strand` set to `TRUE`.
* Default `plotApairOfChrs()` to the first _sequence level_ of the `GBreaks` object.
* `plotApairOfChrs()` now renders regions as blocks and not arrows as the strand
  information is encoded in the color and shape of the comparison bands. It
  also plots in light blue the blocks that are not matched.  Closes issues #3
  and #5.
* Added example data, see ?StructuralVariants and `vignette("Structural Variants",
  package = "GenomicBreaks") for details.

## Bug fixes

* In `flagInversions()`, check if input was sorted.  Closes Issue #2.
* `flagLongShort()` now handles the case where some annotations do not match.

# GenomicBreaks 0.7.1

* Fix namespace and variable bugs that were preventing `reverse()` to run.

# GenomicBreaks 0.7.0

* New `mergeSeqLevels()` function to concatenate sequence features, for instance
  before plotting.
* New `reverse()` method for `GBreaks` and `GRanges` objects.
* Really return a `GBreaks` object in `coalesce_contigs()`.
* Pass extra arguments from `plotApairOfChrs()` to `genoPlotR::plot_gene_map`.

# GenomicBreaks 0.6.0

* Add pairwise chromosome plot functions using the
  [`genoPlotR`](http://genoplotr.r-forge.r-project.org) package,
  see in particular `plotApairOfChrs()`.

# GenomicBreaks 0.5.3

* _pkgdown_ update.

# GenomicBreaks 0.5.2

* Greatly speed up `flagLongShort()`.

# GenomicBreaks 0.5.1

* Correct a variable name bug in `flagLongShort()`.
* Export the `GBreaks()` constructor.

# GenomicBreaks 0.5.0

* Add a `flagLongShort()` function to transfer annotations about chromosome arms.

# GenomicBreaks 0.4.0

* Use alignments of _Saccharomyces paradoxus_ to _S. cerevisiae_ as example data.
* Addition of a `GOC()` function to calculate Gene Order Conservation.

# GenomicBreaks 0.3.0

* Added the `leftInversionGaps()` function.

# GenomicBreaks 0.2.0

* Document the package with [`pkgdown`](https://pkgdown.r-lib.org/).
