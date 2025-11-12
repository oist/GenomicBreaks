# Flip strand of inversions

Flip strand of inversions

## Usage

``` r
flipInversions(gb, detect = TRUE)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- detect:

  Run again
  [`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md)
  if `TRUE`, else reuse the `inv` flag or fail if it is absent.

## Value

Returns the `GBreaks` object in which all ranges that are the central
part of an inversion triplet have their strand orientation flipped. The
`inv` metadata column is then discarded as it is not valid anymore. As
the former inversion triplets are now collinear, new inversions may
possibly found after the object is coalesced again.

## See also

Other Inversion functions:
[`filterDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/filterDoubleInversions.md),
[`filterInversions()`](https://oist.github.io/GenomicBreaks/reference/filterInversions.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`leftInversionGaps()`](https://oist.github.io/GenomicBreaks/reference/leftInversionGaps.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
[`showInversions()`](https://oist.github.io/GenomicBreaks/reference/showInversions.md)

Other modifier functions:
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
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

## Examples

``` r
exampleNestedInversions |> flagInversions()
#> GBreaks object with 5 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query   inv
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle>
#>   [1]     chrA   100-190      + | chrA:100-190 FALSE
#>   [2]     chrA   200-290      - | chrA:400-490  TRUE
#>   [3]     chrA   300-390      + | chrA:300-390 FALSE
#>   [4]     chrA   400-490      - | chrA:200-290 FALSE
#>   [5]     chrA   500-590      + | chrA:500-590 FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
flipInversions(exampleNestedInversions)
#> GBreaks object with 5 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-190      + | chrA:100-190
#>   [2]     chrA   200-290      - | chrA:400-490
#>   [3]     chrA   300-390      - | chrA:300-390
#>   [4]     chrA   400-490      - | chrA:200-290
#>   [5]     chrA   500-590      + | chrA:500-590
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
