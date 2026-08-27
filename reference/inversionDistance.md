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

This function uses several helper functions defined in breakpointGraph.R
(e.g., `extendedPermutation`, `breakpoint_graph`, `hurdles_count`,
`superhurdles_count`, and others) to compute properties of the
breakpoint graph and identify cycles, hurdles and superhurdles. It also
depends on
[`permutationVector()`](https://oist.github.io/GenomicBreaks/reference/permutationVector.md),
which is defined in permutationVector.R.

This algorithm was designed to work in *a single, linear chromosome
alignment*. Although the function still works if the GBreaks object
involves more than one chromosome, the returned value for the minimal
number of inversions will imply in non-usual inversions if different
chromosomes have orthologous regions.

## References

Hannenhalli, Sridhar, and Pavel A. Pevzner. "Transforming cabbage into
turnip: polynomial algorithm for sorting signed permutations by
reversals." Journal of the ACM (JACM) 46.1 (1999): 1-27.

## See also

[`permutationVector`](https://oist.github.io/GenomicBreaks/reference/permutationVector.md)
for generating the permutation vector.

Other Rearrangement distances:
[`breakpointGraphProperties()`](https://oist.github.io/GenomicBreaks/reference/breakpointGraphProperties.md),
[`inversionEstimate_BD()`](https://oist.github.io/GenomicBreaks/reference/inversionEstimate_BD.md)

Other Similarity indexes:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`JC69_distance_allseq()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance_allseq.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`K80_gap_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_gap_distance.md),
[`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md),
[`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md),
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`breakpointGraphProperties()`](https://oist.github.io/GenomicBreaks/reference/breakpointGraphProperties.md),
[`correlation_index()`](https://oist.github.io/GenomicBreaks/reference/correlation_index.md),
[`inversionEstimate_BD()`](https://oist.github.io/GenomicBreaks/reference/inversionEstimate_BD.md),
[`karyotype_index()`](https://oist.github.io/GenomicBreaks/reference/karyotype_index.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md),
[`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md),
[`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md),
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md),
[`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)

## Author

Bruna Fistarol

## Examples

``` r
# Example using a permutation vector directly
# Suppose we want to sort the permutation p = c(1, 3, -2, 4)
inversionDistance(c(1, 3, -2, 4))
#> [1] 2

# Example using a GBreaks object.
# This example is based on Figure 4(a) from Hannehalli and Pevzner (1999). 
# The inversion distance is 8.
inversionDistance(exampleInversionHP1999fig4a)
#> [1] 8
```
