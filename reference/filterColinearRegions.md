# Filter colinear regions

Filter colinear regions

## Usage

``` r
filterColinearRegions(gb, rename = TRUE)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object processed with
  [`flagColinearAlignments`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md).

- rename:

  Replace range names by their numeric order before subsetting.

## Value

Returns the `GBreaks` object in which all ranges that are not the
central part of an inversion triplet have been discarded. If the object
was missing the `colinear` metadata column, return the object after
discarding all of its ranges.

## See also

Other Colinearity functions:
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`chain_contigs()`](https://oist.github.io/GenomicBreaks/reference/chain_contigs.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`dist2next()`](https://oist.github.io/GenomicBreaks/reference/dist2next.md),
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md)

## Examples

``` r
filterColinearRegions(flagColinearAlignments(exampleColinear))
#> GBreaks object with 2 ranges and 2 metadata columns:
#>     seqnames    ranges strand |        query  colinear
#>        <Rle> <IRanges>  <Rle> |    <GRanges> <logical>
#>   1     chrA   100-150      + | chrB:100-150      TRUE
#>   2     chrA   251-300      + | chrB:251-300     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
