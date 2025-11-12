# Coalesce Pairwise Alignments

Coalesce Pairwise Alignments

## Usage

``` r
coalesce_contigs(gb, tol = Inf, minwidth = 0)
```

## Arguments

- gb:

  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object of the pairwise alignment.

- tol:

  Unaligned region of width lesser than or equal to `tol` in both the
  reference and query case will be bridged in coalescing.

- minwidth:

  Remove the intervals whose width smaller than this value.

## Value

Returns a new `GBreaks` object with a reduced number of alignments
fragments due to coalescion. The returned object is sorted ignoring
strand. For convenience during analysis sessions, its *score* is set to
the width of the ranges on the *target* genome.

## Note

Fragmented alignments arising from incorrect basecalls, misassembly or
misalignments can cause us to infer artificial breakpoints

Internally, `coalesce_contigs` uses
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md)
See the examples for details.

## See also

Other Colinearity functions:
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`chain_contigs()`](https://oist.github.io/GenomicBreaks/reference/chain_contigs.md),
[`dist2next()`](https://oist.github.io/GenomicBreaks/reference/dist2next.md),
[`filterColinearRegions()`](https://oist.github.io/GenomicBreaks/reference/filterColinearRegions.md),
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md)

Other modifier functions:
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
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

Charlotte West

Charles Plessy

## Examples

``` r
flagColinearAlignments(exampleColinear3, details = TRUE)
#> GBreaks object with 3 ranges and 10 metadata columns:
#>       seqnames    ranges strand |        query     tfoll     qfoll     tprev
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <integer> <integer> <integer>
#>   [1]     chrA   100-200      + | chrB:100-200         1         1      <NA>
#>   [2]     chrA   201-300      + | chrB:201-300         1         1        -1
#>   [3]     chrA   301-400      + | chrB:301-400      <NA>      <NA>        -1
#>           qprev     t_col     q_col     tdist     qdist  colinear
#>       <integer> <logical> <logical> <numeric> <numeric> <logical>
#>   [1]      <NA>      TRUE      TRUE         1         1      TRUE
#>   [2]        -1      TRUE      TRUE         1         1      TRUE
#>   [3]        -1      <NA>      <NA>       Inf       Inf     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
coalesce_contigs(exampleColinear3)
#> GBreaks object with 1 range and 2 metadata columns:
#>       seqnames    ranges strand |        query     score
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <integer>
#>   [1]     chrA   100-400      + | chrB:100-400       301
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome

# Target range [1] precedes target range [2]
precede(exampleColinear3)
#> [1]  2  3 NA
# Query range [1] precedes query range [2]
precede(exampleColinear3$query)
#> [1]  2  3 NA

# Ranges on the minus strand
gb2 <- exampleColinear3 |> reverse() |> sort(ignore.strand = TRUE)
flagColinearAlignments(gb2, details = TRUE)
#> GBreaks object with 3 ranges and 10 metadata columns:
#>       seqnames    ranges strand |        query     tfoll     qfoll     tprev
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <integer> <integer> <integer>
#>   [1]     chrA   201-300      - | chrB:301-400      <NA>      <NA>         1
#>   [2]     chrA   301-400      - | chrB:201-300        -1        -1         1
#>   [3]     chrA   401-501      - | chrB:100-200        -1        -1      <NA>
#>           qprev     t_col     q_col     tdist     qdist  colinear
#>       <integer> <logical> <logical> <numeric> <numeric> <logical>
#>   [1]         1      TRUE      TRUE         1         1      TRUE
#>   [2]         1      TRUE      TRUE         1         1      TRUE
#>   [3]      <NA>      <NA>      <NA>       Inf       Inf     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
coalesce_contigs(gb2)
#> GBreaks object with 1 range and 2 metadata columns:
#>       seqnames    ranges strand |        query     score
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <integer>
#>   [1]     chrA   201-501      - | chrB:100-400       301
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome

# Target range [1] follows target range [2]
follow(gb2)
#> [1]  2  3 NA
# Or, ignoring strand, target range [1] precedes target range [2]
precede(gb2, ignore.strand = TRUE)
#> [1]  2  3 NA
# Query range [1] follows query range [2]
follow(gb2$query)
#> [1]  2  3 NA

# Coalescing strandless objects
gb3 <- exampleColinear3
gb4 <- gb2
strand(gb4) <- strand(gb3) <- "*"
coalesce_contigs(gb3)
#> GBreaks object with 1 range and 2 metadata columns:
#>       seqnames    ranges strand |        query     score
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <integer>
#>   [1]     chrA   100-400      * | chrB:100-400       301
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
coalesce_contigs(gb4)
#> GBreaks object with 1 range and 2 metadata columns:
#>       seqnames    ranges strand |        query     score
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <integer>
#>   [1]     chrA   201-501      * | chrB:100-400       301
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
