# Construct a `GBreaks` object from a permutation vector

Inverse of
[`permutationVector()`](https://oist.github.io/GenomicBreaks/reference/permutationVector.md):
builds a
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
object of evenly spaced blocks from a signed permutation.

## Usage

``` r
permutationVectorToGBreaks(
  v,
  tsplit = length(v),
  qsplit = length(v),
  width = 90,
  spacing = 100
)
```

## Arguments

- v:

  A signed permutation of `1:n`, as returned by
  [`permutationVector()`](https://oist.github.io/GenomicBreaks/reference/permutationVector.md).

- tsplit, qsplit:

  Number of blocks in each sequence of the *target* and *query* genomes,
  in order. They must sum to `n`.

- width:

  Width of the blocks.

- spacing:

  Distance between the starts of consecutive blocks.

## Value

A `GBreaks` object of length `n`, sorted, such that
[`permutationVector()`](https://oist.github.io/GenomicBreaks/reference/permutationVector.md)
returns `v`.

## Details

Block `i` of the *target* genome is aligned to position `abs(v[i])` of
the *query* genome, on the minus strand if `v[i]` is negative. The
blocks are placed every `spacing` bases on each genome, starting at
position `spacing`, and the genomes are cut into sequences of `tsplit`
and `qsplit` blocks, named `chrT1`, `chrT2`, etc. and `chrQ1`, `chrQ2`,
etc. Sequence lengths are set, so that the object can be passed to
[`reverse()`](https://oist.github.io/GenomicBreaks/reference/reverse.md).

## See also

[`simulateRearrangements()`](https://oist.github.io/GenomicBreaks/reference/simulateRearrangements.md)

## Examples

``` r
# The double inversion of exampleDoubleInversion1
permutationVector(exampleDoubleInversion1)
#> [1]  1  4 -2 -3  5
permutationVectorToGBreaks(c(1, 4, -2, -3, 5))
#> GBreaks object with 5 ranges and 1 metadata column:
#>       seqnames    ranges strand |         query
#>          <Rle> <IRanges>  <Rle> |     <GRanges>
#>   [1]    chrT1   100-189      + | chrQ1:100-189
#>   [2]    chrT1   200-289      + | chrQ1:400-489
#>   [3]    chrT1   300-389      - | chrQ1:200-289
#>   [4]    chrT1   400-489      - | chrQ1:300-389
#>   [5]    chrT1   500-589      + | chrQ1:500-589
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
permutationVectorToGBreaks(c(1, 4, -2, -3, 5)) |> flagDoubleInversions()
#> GBreaks object with 5 ranges and 2 metadata columns:
#>       seqnames    ranges strand |         query       Dbl
#>          <Rle> <IRanges>  <Rle> |     <GRanges> <logical>
#>   [1]    chrT1   100-189      + | chrQ1:100-189      TRUE
#>   [2]    chrT1   200-289      + | chrQ1:400-489     FALSE
#>   [3]    chrT1   300-389      - | chrQ1:200-289     FALSE
#>   [4]    chrT1   400-489      - | chrQ1:300-389     FALSE
#>   [5]    chrT1   500-589      + | chrQ1:500-589     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome

# On two target sequences
permutationVectorToGBreaks(c(1, 4, -2, -3, 5), tsplit = c(2, 3))
#> GBreaks object with 5 ranges and 1 metadata column:
#>       seqnames    ranges strand |         query
#>          <Rle> <IRanges>  <Rle> |     <GRanges>
#>   [1]    chrT1   100-189      + | chrQ1:100-189
#>   [2]    chrT1   200-289      + | chrQ1:400-489
#>   [3]    chrT2   100-189      - | chrQ1:200-289
#>   [4]    chrT2   200-289      - | chrQ1:300-389
#>   [5]    chrT2   300-389      + | chrQ1:500-589
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
```
