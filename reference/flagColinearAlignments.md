# Flag colinear alignments

Flags alignments that are colinear with the next one in sequence order.
The flag is added to the first alignment. Colinearity is defined by the
fact that the next alignment on the *target* genome is the same as the
next alignment on the *query* genome, with strand information being
taken into account.

## Usage

``` r
flagColinearAlignments(gb, tol = Inf, minwidth = 0, details = FALSE)
```

## Arguments

- gb:

  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object of the pairwise alignment.

- tol:

  Unaligned regions larger than this *tolerance* threshold will
  interrupt colinearity.

- minwidth:

  Remove the intervals whose width smaller than this value.

- details:

  Returns more metadata columns if `TRUE`.

## Value

Returns a modified `GBreaks` object with a new `colinear` metadata
column indicating if an alignment is colinear with the next one. If the
`details` option is set to `TRUE`, it will also output the relative
position of the following and previous alignments on the *target* and
*query* genomes (`tfoll`, `tprev`, `qfoll`, `qprev`), a partial
colinearity flag for each genome (`t_col` and `q_col`), and the distance
to the next alignment on each genome (`tdist` and `qdist`).

## Details

Internally, `flagColinearAlignments()` uses the
[`GenomicRanges::precede`](https://rdrr.io/pkg/GenomicRanges/man/nearest-methods.html)
and
[`GenomicRanges::follow`](https://rdrr.io/pkg/GenomicRanges/man/nearest-methods.html)
functions functions to determine what is the *next*. For a given range,
these functions return the index position of the range it precedes or
follows, or `NA` as the first range follows nothing and the last range
precedes nothing. See the examples for details.

## Note

The flags are only valid as long as the Genomic Breaks object is not
modified by removing alignments or sorting them in a different order.

Pay attention that if the `mindwidth` option is passed, some intervals
are discarded from the returned object. This parameter might be removed
in the future if too confusing or useless.

## See also

Other Flagging functions:
[`flagAll()`](https://oist.github.io/GenomicBreaks/reference/flagAll.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagLongShort()`](https://oist.github.io/GenomicBreaks/reference/flagLongShort.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
[`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md)

Other Colinearity functions:
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`chain_contigs()`](https://oist.github.io/GenomicBreaks/reference/chain_contigs.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`dist2next()`](https://oist.github.io/GenomicBreaks/reference/dist2next.md),
[`filterColinearRegions()`](https://oist.github.io/GenomicBreaks/reference/filterColinearRegions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md)

## Author

Charlotte West

Charles Plessy

## Examples

``` r
flagColinearAlignments(exampleColinear3)
#> GBreaks object with 3 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query  colinear
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <logical>
#>   [1]     chrA   100-200      + | chrB:100-200      TRUE
#>   [2]     chrA   201-300      + | chrB:201-300      TRUE
#>   [3]     chrA   301-400      + | chrB:301-400     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
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

# Target range [1] precedes target range [2]
precede(exampleColinear3)
#> [1]  2  3 NA
# Query range [1] precedes query range [2]
precede(exampleColinear3$query)
#> [1]  2  3 NA

# Ranges on the minus strand
gb2 <- exampleColinear3 |> reverse() |> sort(ignore.strand = TRUE)
flagColinearAlignments(gb2)
#> GBreaks object with 3 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query  colinear
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <logical>
#>   [1]     chrA   201-300      - | chrB:301-400      TRUE
#>   [2]     chrA   301-400      - | chrB:201-300      TRUE
#>   [3]     chrA   401-501      - | chrB:100-200     FALSE
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

# Colinearity check on strandless objects
gb3 <- exampleColinear3
gb4 <- gb2
strand(gb4) <- strand(gb3) <- "*"
flagColinearAlignments(gb3)
#> GBreaks object with 3 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query  colinear
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <logical>
#>   [1]     chrA   100-200      * | chrB:100-200      TRUE
#>   [2]     chrA   201-300      * | chrB:201-300      TRUE
#>   [3]     chrA   301-400      * | chrB:301-400     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
flagColinearAlignments(gb4)
#> GBreaks object with 3 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query  colinear
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <logical>
#>   [1]     chrA   201-300      * | chrB:301-400      TRUE
#>   [2]     chrA   301-400      * | chrB:201-300      TRUE
#>   [3]     chrA   401-501      * | chrB:100-200     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome

# Ranges that should not coalesce because they are not
# ordered properly
flagColinearAlignments(exampleNotColinear)
#> GBreaks object with 2 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query  colinear
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <logical>
#>   [1]     chrA   100-150      + | chrB:201-251     FALSE
#>   [2]     chrA   251-300      + |  chrB:50-100     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
