# Bridge regions

Maps unaligned regions of the *target* to the *query* genome when they
are flanked by colinear regions.

## Usage

``` r
bridgeRegions(gb)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

## Value

Returns a new `GBreaks` object of shorter length. Its length will be
zero if no bridge regions were found.

## Note

Because some aligned regions can be directly adjacent (no gaps), the
returned `GBreaks` object may contain ranges of width zero, where the
*start* coordinate is 1 nucleotide higher than the *end* coordinate.

## References

Bridge regions have also been called “simultaneous gaps” in the
comparison of the mouse and human genomes by Kent WJ, Baertsch R,
Hinrichs A, Miller W, Haussler D. (*Evolution's cauldron: duplication,
deletion, and rearrangement in the mouse and human genomes.* Proc Natl
Acad Sci U S A. 2003;100(20):11484-11489. doi:10.1073/pnas.1932072100)

## See also

Other Colinearity functions:
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`chain_contigs()`](https://oist.github.io/GenomicBreaks/reference/chain_contigs.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`dist2next()`](https://oist.github.io/GenomicBreaks/reference/dist2next.md),
[`filterColinearRegions()`](https://oist.github.io/GenomicBreaks/reference/filterColinearRegions.md),
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md)

Other modifier functions:
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`forceSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/forceSeqLengths.md),
[`guessSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/guessSeqLengths.md),
[`keepLongestPair()`](https://oist.github.io/GenomicBreaks/reference/keepLongestPair.md),
[`matchPairs()`](https://oist.github.io/GenomicBreaks/reference/matchPairs.md),
[`mergeSeqLevels()`](https://oist.github.io/GenomicBreaks/reference/mergeSeqLevels.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
[`removeTranslocations()`](https://oist.github.io/GenomicBreaks/reference/removeTranslocations.md),
[`reverse()`](https://oist.github.io/GenomicBreaks/reference/reverse.md),
[`splitSeqLevel()`](https://oist.github.io/GenomicBreaks/reference/splitSeqLevel.md),
[`swap()`](https://oist.github.io/GenomicBreaks/reference/swap.md)

## Author

Charles Plessy

## Examples

``` r
exampleColinear5
#> GBreaks object with 5 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-190      + | chrB:100-190
#>   [2]     chrA   200-290      + | chrB:200-290
#>   [3]     chrA   300-390      + | chrB:300-390
#>   [4]     chrA   400-490      + | chrB:400-490
#>   [5]     chrA   500-590      + | chrB:500-590
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
bridgeRegions(exampleColinear5)
#> GBreaks object with 4 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   191-199      * | chrB:191-199
#>   [2]     chrA   291-299      * | chrB:291-299
#>   [3]     chrA   391-399      * | chrB:391-399
#>   [4]     chrA   491-499      * | chrB:491-499
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome

# Note the zero-width ranges when aligned regions are directly adjacent.
exampleColinear3
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-200      + | chrB:100-200
#>   [2]     chrA   201-300      + | chrB:201-300
#>   [3]     chrA   301-400      + | chrB:301-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
bridgeRegions(exampleColinear3)
#> GBreaks object with 2 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   201-200      * | chrB:201-200
#>   [2]     chrA   301-300      * | chrB:301-300
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
