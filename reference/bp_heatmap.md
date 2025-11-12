# Breakpoint Associated Heatmaps

This function generates a heatmap of the specified pattern, over
breakpoints aligned at the centre of the plot. The result can be piped
into
[`heatmaps::smoothHeatmap`](https://rdrr.io/pkg/heatmaps/man/smoothHeatmap.html)
and then
[`heatmaps::plotHeatmapList`](https://rdrr.io/pkg/heatmaps/man/plotHeatmapList.html)
or
[`heatmaps::plotHeatmapMeta`](https://rdrr.io/pkg/heatmaps/man/plotHeatmapMeta.html)

## Usage

``` r
bp_heatmap(gr, window, pattern, ...)
```

## Arguments

- gr:

  `GRanges` object containing pairwise alignment

- window:

  Range over which heatmap characteristic is plotted. Breakpoints will
  be aligned at the center of this.

- pattern:

  Character string of desired pattern/characteristic to be plotted on
  heatmap.

- ...:

  Pass other arguments to
  [`get_bps`](https://oist.github.io/GenomicBreaks/reference/get_bps.md).

## Value

Returns a
[`heatmaps::Heatmap`](https://rdrr.io/pkg/heatmaps/man/Heatmap.html)
object of `pattern` around centred breakpoints.

## Note

The `GRanges` object is expected to have a *sequence information* (see
[`Seqinfo::seqinfo`](https://rdrr.io/pkg/Seqinfo/man/seqinfo.html)) that
allows the retrieval its corresponding `BSgenome` object via the
[`BSgenome::getBSgenome`](https://rdrr.io/pkg/BSgenome/man/available.genomes.html)
function.

## See also

Other plot functions:
[`bp_pair_analysis()`](https://oist.github.io/GenomicBreaks/reference/bp_pair_analysis.md),
[`feature_coverage()`](https://oist.github.io/GenomicBreaks/reference/feature_coverage.md),
[`makeOxfordPlots()`](https://oist.github.io/GenomicBreaks/reference/makeOxfordPlots.md),
[`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)

Other heatmap functions:
[`bp_pair_analysis()`](https://oist.github.io/GenomicBreaks/reference/bp_pair_analysis.md),
[`feature_coverage()`](https://oist.github.io/GenomicBreaks/reference/feature_coverage.md)

## Author

Charlotte West

Charles Plessy

## Examples

``` r
# The plot makes no sense, but that is the best example I have at the moment.
exdata_Sac <- system.file("extdata/SacCer3__SacPar.gff3.gz", package = "GenomicBreaks")
gb <- load_genomic_breaks(exdata_Sac, BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae)
bp_heatmap(gb, 200, 'GC', dir = "left") |>
  heatmaps::smoothHeatmap() |> heatmaps::plotHeatmapList()
#> 
#> Attaching package: ‘Biostrings’
#> The following object is masked from ‘package:GenomicBreaks’:
#> 
#>     pairwiseAlignment
#> The following object is masked from ‘package:base’:
#> 
#>     strsplit
#> 
#> Calculating kernel density...
#> plotting heatmap GC

```
