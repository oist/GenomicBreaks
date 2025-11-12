# Compute all flags

Flag ranges that are at the beginning of a colinear duplet, or an
inversion or a translocation triplet.

## Usage

``` r
flagAll(gb)
```

## Arguments

- gb:

  A `GBreaks` object.

## Value

Returns the `GBreaks` object with an extra `flag` metadata column.

## See also

Other Flagging functions:
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagLongShort()`](https://oist.github.io/GenomicBreaks/reference/flagLongShort.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md)

## Examples

``` r
flagAll(exampleInversion)
#> GBreaks object with 3 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query        flag
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <character>
#>   [1]     chrA   100-190      + | chrB:100-190         Inv
#>   [2]     chrA   210-291      - | chrB:210-291        <NA>
#>   [3]     chrA   301-400      + | chrB:301-400        <NA>
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
