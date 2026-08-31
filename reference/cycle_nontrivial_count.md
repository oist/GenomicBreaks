# Number of Non-Trivial Cycles

This function computes the number of **non-trivial** cycles in a
breakpoint graph. Non-trivial cycles are cycles with more than 2
vertices, whereas trivial cycles are cycles with 2 vertices that
represent common adjacencies.

## Usage

``` r
cycle_nontrivial_count(g)
```

## Arguments

- g:

  The breakpoint graph.

## Value

The number of cycles in a breakpoint graph.

## See also

Other Breakpoint graph functions:
[`bp_count()`](https://oist.github.io/GenomicBreaks/reference/bp_count.md),
[`breakpointGraphProperties()`](https://oist.github.io/GenomicBreaks/reference/breakpointGraphProperties.md),
[`breakpoint_graph()`](https://oist.github.io/GenomicBreaks/reference/breakpoint_graph.md),
[`components_graph()`](https://oist.github.io/GenomicBreaks/reference/components_graph.md),
[`cycle_count()`](https://oist.github.io/GenomicBreaks/reference/cycle_count.md),
[`extendedPermutation()`](https://oist.github.io/GenomicBreaks/reference/extendedPermutation.md),
[`hurdles_count()`](https://oist.github.io/GenomicBreaks/reference/hurdles_count.md),
[`is_fortress()`](https://oist.github.io/GenomicBreaks/reference/is_fortress.md),
[`is_interleaving()`](https://oist.github.io/GenomicBreaks/reference/is_interleaving.md),
[`superhurdles_count()`](https://oist.github.io/GenomicBreaks/reference/superhurdles_count.md)

## Author

Bruna Fistarol
