# Chain contigs in pairwise alignments

Description TBW

## Usage

``` r
chain_contigs(gb, tol = Inf, drop.within = FALSE)
```

## Arguments

- gb:

  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object of the pairwise alignment.

- tol:

  Unaligned region of width lesser than or equal to `tol` in both the
  reference and query case will be bridged in coalescing.

- drop.within:

  Discard the chains that are included in other chains.

## Value

TBW

## See also

Other Colinearity functions:
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`dist2next()`](https://oist.github.io/GenomicBreaks/reference/dist2next.md),
[`filterColinearRegions()`](https://oist.github.io/GenomicBreaks/reference/filterColinearRegions.md),
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md)

## Author

Charles Plessy

## Examples

``` r
chain_contigs(exampleInversion)
#> GBreaks object with 2 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query     score
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <integer>
#>   [1]     chrA   100-400      + | chrB:100-400       301
#>   [2]     chrA   210-291      - | chrB:210-291        82
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
