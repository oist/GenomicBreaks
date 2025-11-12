# Flag long and short arms

*Oikopleura* genomes are special as the long and short arms of their
chromosomes have different properties such as `GC` or repeat content. It
can be useful to know if genomic region belongs to long or a short arm.

## Usage

``` r
flagLongShort(gr, annot, select = c("first", "last", "arbitrary"))
```

## Arguments

- gr:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  or a simple
  [`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object

- annot:

  A `GRanges` file containing the coordinate of arms and their nature
  (such as `long`, `short`, `XSR` or `YSR`) in a `Type` *metadata
  column*.

- select:

  One method among `first`, `last` and `arbitrary`, to decide on how to
  handle the regions that match both arms.

## Value

Returns a modified version of the object with an extra *metadata
column*, `Arm` in which the information from the annotation file's
`Type` column was transferred. See the examples below and in the manual
of
[`IRanges::findOverlaps`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)
for details. on how regions that match both arms are handled.

## See also

Other Flagging functions:
[`flagAll()`](https://oist.github.io/GenomicBreaks/reference/flagAll.md),
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md)

## Author

Charles Plessy

## Examples

``` r
annot <- GRanges(c("chrA:1-140", "chrA:150-500"), Type = c("short", "long"))

flagLongShort(exampleColinear, annot)
#> GBreaks object with 2 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query      Arm
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <factor>
#>   [1]     chrA   100-150      + | chrB:100-150    short
#>   [2]     chrA   251-300      + | chrB:251-300    long 
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
flagLongShort(exampleColinear, annot, select = "last")
#> GBreaks object with 2 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query      Arm
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <factor>
#>   [1]     chrA   100-150      + | chrB:100-150     long
#>   [2]     chrA   251-300      + | chrB:251-300     long
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
flagLongShort(exampleColinear, annot, select = "arbitrary")
#> GBreaks object with 2 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query      Arm
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <factor>
#>   [1]     chrA   100-150      + | chrB:100-150    short
#>   [2]     chrA   251-300      + | chrB:251-300    long 
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
