# Remove translocations

Remove translocations

## Usage

``` r
removeTranslocations(gb, detect = TRUE)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- detect:

  Run again
  [`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md)
  if `TRUE`, else reuse the `tra` and `trs` flags. If they are absent,
  an empty object is returned.

## Value

Returns the `GBreaks` object after removing the moved ranges (see
[`filterTranslocations()`](https://oist.github.io/GenomicBreaks/reference/filterTranslocations.md))
and coalescing the resulting collinear ranges.

## See also

Other Translocation functions:
[`filterTranslocations()`](https://oist.github.io/GenomicBreaks/reference/filterTranslocations.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
[`showTranslocations()`](https://oist.github.io/GenomicBreaks/reference/showTranslocations.md)

Other modifier functions:
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`forceSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/forceSeqLengths.md),
[`guessSeqLengths()`](https://oist.github.io/GenomicBreaks/reference/guessSeqLengths.md),
[`keepLongestPair()`](https://oist.github.io/GenomicBreaks/reference/keepLongestPair.md),
[`matchPairs()`](https://oist.github.io/GenomicBreaks/reference/matchPairs.md),
[`mergeSeqLevels()`](https://oist.github.io/GenomicBreaks/reference/mergeSeqLevels.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
[`reverse()`](https://oist.github.io/GenomicBreaks/reference/reverse.md),
[`splitSeqLevel()`](https://oist.github.io/GenomicBreaks/reference/splitSeqLevel.md),
[`swap()`](https://oist.github.io/GenomicBreaks/reference/swap.md)

## Examples

``` r
removeTranslocations(exampleTranslocation)
#> GBreaks object with 1 range and 2 metadata columns:
#>       seqnames    ranges strand |        query     score
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <integer>
#>   [1]     chrA   100-400      + | chrB:100-400       301
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
removeTranslocations(swap(exampleTranslocation, sort = TRUE))
#> GBreaks object with 1 range and 2 metadata columns:
#>       seqnames    ranges strand |        query     score
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <integer>
#>   [1]     chrB   100-400      + | chrA:100-400       301
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
```
