# Inversion Distance

Computes the minimal number of inversions required to sort a signed
permutation using the Hannenhalli and Pevzner algorithm.

## Usage

``` r
inversionDistance(x)
```

## Arguments

- x:

  Either a GBreaks object or a signed permutation vector. If `x` is a
  GBreaks object, a permutation vector will be extracted using
  [`permutationVector`](https://oist.github.io/GenomicBreaks/reference/permutationVector.md).

## Value

An integer: the minimal number of inversions needed to sort the
permutation.

## Details

This function uses several internal helper functions (e.g.,
`extendedPermutation`, `breakpoint_graph`, `hurdles_count`,
`superhurdles_count`, and others) to compute properties of the
breakpoint graph and identify cycles, hurdles and superhurdles. It also
depends on
[`permutationVector()`](https://oist.github.io/GenomicBreaks/reference/permutationVector.md),
which is defined in another .R file.

This algorithm was designed to work in *a single, linear chromosome
alignment*. Although the function still works if the GBreaks object
involves more than one chromosome, the returned value for the minimal
number of inversions will imply in non-usual inversions if different
chromosomes have orthologous regions.

## See also

[`permutationVector`](https://oist.github.io/GenomicBreaks/reference/permutationVector.md)
for generating the permutation vector.

## Author

Bruna Fistarol

## Examples

``` r
if (FALSE) { # \dontrun{
# Example using a permutation vector directly
# Suppose we want to sort the permutation p = c(1, 3, -2, 4)
inversionDistance(c(1, 3, -2, 4))

# Example using a GBreaks object
gr         <- GRanges(c("A:10-19", "A:20-29", "A:30-39", "A:40-49", "A:50-59", "A:60-69"))
strand(gr) <- c(            "+"  ,     "-"  ,     "-"  ,     "-"  ,     "+"  ,     "+"   )
gr$query   <- GRanges(c("A:10-19", "A:30-39", "A:20-29", "A:40-49", "A:60-69", "A:50-59"))
gb         <- GBreaks(gr)
inversionDistance(gb)} # }
```
