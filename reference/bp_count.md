# Number of Breakpoints

This function computes the number of breakpoints given an extended
permutation.

## Usage

``` r
bp_count(p_extended)
```

## Arguments

- p_extended:

  The extended permutation.

## Value

The number of breakpoints between two genomes.

## See also

Other Breakpoint graph functions:
[`breakpointGraphProperties()`](https://oist.github.io/GenomicBreaks/reference/breakpointGraphProperties.md),
[`breakpoint_graph()`](https://oist.github.io/GenomicBreaks/reference/breakpoint_graph.md),
[`components_graph()`](https://oist.github.io/GenomicBreaks/reference/components_graph.md),
[`cycle_count()`](https://oist.github.io/GenomicBreaks/reference/cycle_count.md),
[`cycle_nontrivial_count()`](https://oist.github.io/GenomicBreaks/reference/cycle_nontrivial_count.md),
[`extendedPermutation()`](https://oist.github.io/GenomicBreaks/reference/extendedPermutation.md),
[`hurdles_count()`](https://oist.github.io/GenomicBreaks/reference/hurdles_count.md),
[`is_fortress()`](https://oist.github.io/GenomicBreaks/reference/is_fortress.md),
[`is_interleaving()`](https://oist.github.io/GenomicBreaks/reference/is_interleaving.md),
[`superhurdles_count()`](https://oist.github.io/GenomicBreaks/reference/superhurdles_count.md)

## Author

Bruna Fistarol

## Examples

``` r
permutationVector(exampleInversion) |>
  GenomicBreaks:::extendedPermutation() |>
  GenomicBreaks:::bp_count()
#> [1] 2
```
