# Ordering permutation for *query* relative to *target*

Computes a permutation which rearranges the sequence levels of the
*query* genome so that it reflects the order in which we see the matches
on a pairwise plot.

## Usage

``` r
orderQuerySeqLevels(gb)

orderQuerySeqLevels_DF_GR(DF, gr)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object that has only one sequence level in use on the *target* genome.

- DF:

  A
  [`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  object representing the *target* genome.

- gr:

  A
  [`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object representing the *query* genome.

## Value

Returns an integer vector of order permutations for the sequence levels
of the *query* genome.

## Details

This is done by computing the average midpoint position of the *query*
ranges on the *target* genome for each seqlevel, matched by the width of
the *target* ranges, so that long matches have more importance. This
avoids spurious ordering due to short matches in the subtelomeric
regions.

## Note

The order only makes sense relative to a single sequence level of the
*target* genome, so the function will stop with error if there was more
than one.

## Author

Charles Plessy

## Examples

``` r
gb       <- GRanges(c("chr1:101-180:+", "chr1:201-300:+",  "chr1:320-400:+"))
gb$query <- GRanges(c( "cgt8:1-100",      "ctg3:1-100",    "ctg5:1-100"))
gb <- GBreaks(gb)
# Sort alphabetically as if they were inherited from a BSgenome object
seqlevels(gb$query) <- c("ctg3", "ctg5", "cgt8")
seqlevels(gb$query)
#> [1] "ctg3" "ctg5" "cgt8"
# Sort by match positions on the target genome.
seqlevels(gb$query) <- seqlevels(gb$query)[orderQuerySeqLevels(gb)]
seqlevels(gb$query)
#> [1] "cgt8" "ctg3" "ctg5"
```
