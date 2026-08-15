# Extended Permutation

This function transforms a permutation with `N` elements into a
permutation with `2N` elements. The idea is that each element in the
input permutation corresponds to an aligned block (or gene, conserved
region, etc.). This function creates two elements for each original
element: one representing the start of the block and the other
representing the end.

## Usage

``` r
extendedPermutation(p)
```

## Arguments

- p:

  A signed permutation vector.

## Value

The extended permutation.

## See also

Other Breakpoint graph functions:
[`bp_count()`](https://oist.github.io/GenomicBreaks/reference/bp_count.md),
[`breakpointGraphProperties()`](https://oist.github.io/GenomicBreaks/reference/breakpointGraphProperties.md),
[`breakpoint_graph()`](https://oist.github.io/GenomicBreaks/reference/breakpoint_graph.md),
[`components_graph()`](https://oist.github.io/GenomicBreaks/reference/components_graph.md),
[`cycle_count()`](https://oist.github.io/GenomicBreaks/reference/cycle_count.md),
[`hurdles_count()`](https://oist.github.io/GenomicBreaks/reference/hurdles_count.md),
[`is_fortress()`](https://oist.github.io/GenomicBreaks/reference/is_fortress.md),
[`is_interleaving()`](https://oist.github.io/GenomicBreaks/reference/is_interleaving.md),
[`superhurdles_count()`](https://oist.github.io/GenomicBreaks/reference/superhurdles_count.md)

## Author

Bruna Fistarol
