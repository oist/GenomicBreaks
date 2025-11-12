# Make Oxford Plots

Takes a
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
object and prepares an Oxford (macrosynteny) of the coordinates of the
*query* ranges against the *target* ranges after concatenating them.

## Usage

``` r
makeOxfordPlots(
  gb,
  sp1Name = "target",
  sp2Name = "query",
  sp1ChrArms = NULL,
  sp2ChrArms = NULL,
  type = c("line", "point", "none"),
  size = 1,
  diag = TRUE,
  col = c("seqnames", "strand", "score")
)
```

## Arguments

- gb:

  A `GBreaks` object

- sp1Name:

  Name of the first species (default: sp1)

- sp2Name:

  Name of the second species (default: sp2)

- sp1ChrArms:

  A `GBreaks` object of chromosome arms in sp1 genome

- sp2ChrArms:

  A `GBreaks` object of chromosome arms in sp2 genome

- type:

  The type of the plot (`point`, `line` or `none`)

- size:

  The size of the plotted dots or segments.

- diag:

  Diagonalise the plot by reordering *query* sequence levels.

- col:

  Colour of the plotted dots or lines by `seqnames`, `strand` or
  `score`.

## Value

Returns a `ggplot2` object that can be further modified using the `+`
operator. Use `type = 'none'` to receive an object without *geom* layer.

## See also

Other plot functions:
[`bp_heatmap()`](https://oist.github.io/GenomicBreaks/reference/bp_heatmap.md),
[`bp_pair_analysis()`](https://oist.github.io/GenomicBreaks/reference/bp_pair_analysis.md),
[`feature_coverage()`](https://oist.github.io/GenomicBreaks/reference/feature_coverage.md),
[`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)

## Author

Aleksandra Bliznina

Charles Plessy

## Examples

``` r
makeOxfordPlots(exampleTranslocation)

makeOxfordPlots(exampleTranslocation, type = "p")

```
