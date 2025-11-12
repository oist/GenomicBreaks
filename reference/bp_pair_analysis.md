# Breakpoint Pairwise Comparison using Reference Genome as Coordinate System

This function takes two
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
objects, with the same *target* genome. A heatmap is produced, where the
breakpoints of *query* 1 are the centre point reference, and the
breakpoints of *query* 2 that fall within the window are projected onto
it.

## Usage

``` r
bp_pair_analysis(gb1, gb2, window, label)
```

## Arguments

- gb1:

  `GBreaks` object of the alignment between the *target* genome and
  *query* genome 1.

- gb2:

  `GBreaks` object of the alignment between the *target* genome and
  *query* genome 2.

- window:

  range around query1 breakpoints of which query2 breakpoints are
  projected on to. Should be even number

- label:

  string which will be the label on the graph

## Value

Returns a
[`heatmaps::Heatmap`](https://rdrr.io/pkg/heatmaps/man/Heatmap.html)
object of `pattern` that can be piped into
[`heatmaps::smoothHeatmap`](https://rdrr.io/pkg/heatmaps/man/smoothHeatmap.html)
and then
[`heatmaps::plotHeatmapList`](https://rdrr.io/pkg/heatmaps/man/plotHeatmapList.html)
or
[`heatmaps::plotHeatmapMeta`](https://rdrr.io/pkg/heatmaps/man/plotHeatmapMeta.html).

## See also

Other plot functions:
[`bp_heatmap()`](https://oist.github.io/GenomicBreaks/reference/bp_heatmap.md),
[`feature_coverage()`](https://oist.github.io/GenomicBreaks/reference/feature_coverage.md),
[`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md),
[`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)

Other heatmap functions:
[`bp_heatmap()`](https://oist.github.io/GenomicBreaks/reference/bp_heatmap.md),
[`feature_coverage()`](https://oist.github.io/GenomicBreaks/reference/feature_coverage.md)

## Author

Charlotte West
