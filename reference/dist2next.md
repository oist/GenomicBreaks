# Distances to next ranges

Calculates the distance to the next range in the *target* and *query*
genome.

## Usage

``` r
dist2next(x, step = 1, ignore.strand = FALSE)

# S4 method for class 'GRanges'
dist2next(x, step = 1, ignore.strand = FALSE)

# S4 method for class 'GBreaks'
dist2next(x, step = 1, ignore.strand = FALSE)
```

## Arguments

- x:

  A
  [`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  or a
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- step:

  Distance to the n^(th) block (default: first).

- ignore.strand:

  Calculate distance for ranges on different strands.

## Value

For `GRanges`, returns the object with one extra metadata colums,
`dist`, and for `GBreaks` two extra columns `tdist` and `qdist`,
containing the distance to the next range in the reference and query
genomes respectively.

## Details

The distance is defined by the
[`GenomicRanges::distance`](https://rdrr.io/pkg/GenomicRanges/man/nearest-methods.html)
function. Note that because the last range has not next neighbor, the
last value is set to `Inf` arbitrarily.

## See also

Other Colinearity functions:
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`chain_contigs()`](https://oist.github.io/GenomicBreaks/reference/chain_contigs.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`filterColinearRegions()`](https://oist.github.io/GenomicBreaks/reference/filterColinearRegions.md),
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md)

## Examples

``` r
dist2next(exampleInversion)
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query     tdist     qdist
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <numeric> <numeric>
#>   [1]     chrA   100-190      + | chrB:100-190        NA        20
#>   [2]     chrA   210-291      - | chrB:210-291        NA        10
#>   [3]     chrA   301-400      + | chrB:301-400       Inf       Inf
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
dist2next(granges(exampleInversion))
#> GRanges object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |      dist
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]     chrA   100-190      + |        NA
#>   [2]     chrA   210-291      - |        NA
#>   [3]     chrA   301-400      + |       Inf
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
dist2next(exampleInversion, ignore.strand = TRUE)
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query     tdist     qdist
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <numeric> <numeric>
#>   [1]     chrA   100-190      + | chrB:100-190        20        20
#>   [2]     chrA   210-291      - | chrB:210-291        10        10
#>   [3]     chrA   301-400      + | chrB:301-400       Inf       Inf
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
dist2next(exampleInversion - 20, ignore.strand = TRUE)
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query     tdist     qdist
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <numeric> <numeric>
#>   [1]     chrA   120-170      + | chrB:100-190        60        20
#>   [2]     chrA   230-271      - | chrB:210-291        50        10
#>   [3]     chrA   321-380      + | chrB:301-400       Inf       Inf
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
dist2next(exampleInversion, 2)
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query     tdist     qdist
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <numeric> <numeric>
#>   [1]     chrA   100-190      + | chrB:100-190       111       111
#>   [2]     chrA   210-291      - | chrB:210-291       Inf       Inf
#>   [3]     chrA   301-400      + | chrB:301-400       Inf       Inf
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
