# Breakpoint Graph

This function creates a breakpoint graph given an extended permutation.

## Usage

``` r
breakpoint_graph(p_extended)
```

## Arguments

- p_extended:

  The extended permutation.

## Value

The breakpoint graph.

## References

Hannenhalli, Sridhar, and Pavel A. Pevzner. "Transforming cabbage into
turnip: polynomial algorithm for sorting signed permutations by
reversals." Journal of the ACM (JACM) 46.1 (1999): 1-27.

## See also

Other Breakpoint graph functions:
[`bp_count()`](https://oist.github.io/GenomicBreaks/reference/bp_count.md),
[`breakpointGraphProperties()`](https://oist.github.io/GenomicBreaks/reference/breakpointGraphProperties.md),
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
  GenomicBreaks:::breakpoint_graph()
#> IGRAPH cf88c0b U--- 8 4 -- 
#> + attr: color (e/c), unoriented (e/n)
#> + edges from cf88c0b:
#> [1] 3--5 3--4 4--6 5--6
```
