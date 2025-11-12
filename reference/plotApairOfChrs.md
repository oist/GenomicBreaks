# Plots a pair of chromosomes from a `GBreaks` object.

Plots the mapping between a pair of chromosomes using the
[`genoPlotR::plot_gene_map`](https://rdrr.io/pkg/genoPlotR/man/plot_gene_map.html)
function. One sequence feature (usually a chromosome) is selected from
the *target* genome. The corresponding sequence feature on the *query*
genome is either chosen automatically (being the one with the largest
fraction of the mappings) or given as a parameter.

## Usage

``` r
plotApairOfChrs(
  gb,
  chrT = NULL,
  chrQ = NULL,
  dna_seg_scale = TRUE,
  dna_seg_labels = NULL,
  ...
)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- chrT:

  A sequence name on the *target* genome. Defaults to the first
  *sequence level* of the `gb` object.

- chrQ:

  (Optional) A sequence name on the *query* genome. Defaults to the
  longest cumulative match on `chrT`.

- dna_seg_scale:

  Plot coordinats and a scale bar (defaults to `TRUE`).

- dna_seg_labels:

  A character vector to override default labels for sequence names.

- ...:

  Further arguments are passed to `plot_gene_map`.

## Value

Plots to the active device and returns and `oma::oma_layout` object like
the
[`genoPlotR::plot_gene_map`](https://rdrr.io/pkg/genoPlotR/man/plot_gene_map.html)
function.

## Note

In this kind of plot, the ranges from the selected sequence on the
*target* genome that are not mapped to the *query* sequence are not
displayed.

## See also

Other genoPlotR functions:
[`gb2comp()`](https://oist.github.io/GenomicBreaks/reference/gb2comp.md),
[`gb2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gb2dna_seg.md),
[`gb2xlim()`](https://oist.github.io/GenomicBreaks/reference/gb2xlim.md),
[`gr2dna_seg()`](https://oist.github.io/GenomicBreaks/reference/gr2dna_seg.md)

Other plot functions:
[`bp_heatmap()`](https://oist.github.io/GenomicBreaks/reference/bp_heatmap.md),
[`bp_pair_analysis()`](https://oist.github.io/GenomicBreaks/reference/bp_pair_analysis.md),
[`feature_coverage()`](https://oist.github.io/GenomicBreaks/reference/feature_coverage.md),
[`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md)

Other Structural variants:
[`StructuralVariants`](https://oist.github.io/GenomicBreaks/reference/StructuralVariants.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md)

## Author

Charles Plessy

## Examples

``` r
plotApairOfChrs(exampleInversion)

plotApairOfChrs(exampleDeletion)


# Labels can be overriden.
plotApairOfChrs(exampleDeletion, dna_seg_labels = c("over...", "...ridden"))

```
