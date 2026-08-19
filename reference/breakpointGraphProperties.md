# Breakpoint graph properties

Computes the breakpoint graph of a pair of chromosomes and returns key
properties of the graph, which are used in various rearrangement
estimates.

## Usage

``` r
breakpointGraphProperties(gb)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

## Value

A list containing the following properties of the breakpoint graph: 1.
`N` : the total number of places where a breakpoint could occur (which
is the same as the `number_aligned_blocks + 1`); 2. `nbBreakpoints` :
the total number of breakpoints in the extended permutation; 3.
`nbCycles` : the total number of cycles in the breakpoint graph.

## References

Hannenhalli, Sridhar, and Pavel A. Pevzner. "Transforming cabbage into
turnip: polynomial algorithm for sorting signed permutations by
reversals." Journal of the ACM (JACM) 46.1 (1999): 1-27.

## See also

Other Breakpoint graph functions:
[`bp_count()`](https://oist.github.io/GenomicBreaks/reference/bp_count.md),
[`breakpoint_graph()`](https://oist.github.io/GenomicBreaks/reference/breakpoint_graph.md),
[`components_graph()`](https://oist.github.io/GenomicBreaks/reference/components_graph.md),
[`cycle_count()`](https://oist.github.io/GenomicBreaks/reference/cycle_count.md),
[`extendedPermutation()`](https://oist.github.io/GenomicBreaks/reference/extendedPermutation.md),
[`hurdles_count()`](https://oist.github.io/GenomicBreaks/reference/hurdles_count.md),
[`is_fortress()`](https://oist.github.io/GenomicBreaks/reference/is_fortress.md),
[`is_interleaving()`](https://oist.github.io/GenomicBreaks/reference/is_interleaving.md),
[`superhurdles_count()`](https://oist.github.io/GenomicBreaks/reference/superhurdles_count.md)

Other Rearrangement distances:
[`inversionDistance()`](https://oist.github.io/GenomicBreaks/reference/inversionDistance.md),
[`inversionEstimate_BD()`](https://oist.github.io/GenomicBreaks/reference/inversionEstimate_BD.md)

Other Similarity indexes:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`K80_gap_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_gap_distance.md),
[`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md),
[`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md),
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`correlation_index()`](https://oist.github.io/GenomicBreaks/reference/correlation_index.md),
[`inversionDistance()`](https://oist.github.io/GenomicBreaks/reference/inversionDistance.md),
[`inversionEstimate_BD()`](https://oist.github.io/GenomicBreaks/reference/inversionEstimate_BD.md),
[`karyotype_index()`](https://oist.github.io/GenomicBreaks/reference/karyotype_index.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md),
[`slidingWindow()`](https://oist.github.io/GenomicBreaks/reference/slidingWindow.md),
[`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md),
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md),
[`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)

## Author

Priscila Biller
