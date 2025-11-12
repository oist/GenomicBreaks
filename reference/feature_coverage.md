# Feature coverage heatmap

Wraps the
[`heatmaps::CoverageHeatmap`](https://rdrr.io/pkg/heatmaps/man/CoverageHeatmap.html)
function to produce a heatmap centred on the boundaries of *genomic
ranges*. Assuming that these ranges are a
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
object, then the boundaries approximate breakpoints.

## Usage

``` r
feature_coverage(gr, feat, window, label, ...)
```

## Arguments

- gr:

  A
  [`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object.

- feat:

  A `GRanges` object representing the feature of interest.

- window:

  Window over which to observe feature coverage.

- label:

  Label for the plot.

- ...:

  Other arguments passed to
  [`get_bps`](https://oist.github.io/GenomicBreaks/reference/get_bps.md)
  in order to select the boundaries, their order and their orientation.

## Value

Returns a
[`heatmaps::Heatmap`](https://rdrr.io/pkg/heatmaps/man/Heatmap.html)
object that can be piped into
[`heatmaps::smoothHeatmap`](https://rdrr.io/pkg/heatmaps/man/smoothHeatmap.html)
and then
[`heatmaps::plotHeatmapList`](https://rdrr.io/pkg/heatmaps/man/plotHeatmapList.html)
or
[`heatmaps::plotHeatmapMeta`](https://rdrr.io/pkg/heatmaps/man/plotHeatmapMeta.html).

## See also

Other plot functions:
[`bp_heatmap()`](https://oist.github.io/GenomicBreaks/reference/bp_heatmap.md),
[`bp_pair_analysis()`](https://oist.github.io/GenomicBreaks/reference/bp_pair_analysis.md),
[`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md),
[`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)

Other heatmap functions:
[`bp_heatmap()`](https://oist.github.io/GenomicBreaks/reference/bp_heatmap.md),
[`bp_pair_analysis()`](https://oist.github.io/GenomicBreaks/reference/bp_pair_analysis.md)
