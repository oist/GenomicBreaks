# Flag twin inversions

Flag ranges that start a quadruplet that would be colinear if the two
central pairs were aligned to the opposite strand. This implies
breakpoint reuse.

## Usage

``` r
flagTwinInversions(gb, details = FALSE)
```

## Arguments

- gb:

  A
  [`GBreaks()`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- details:

  Returns more metadata columns if `TRUE`.

## Value

Returns the `GBreaks` object with an extra `twi` metadata column.

## Details

Here is a trivial example of a twin inversion.

    ┌──────────────┬──────────────┬──────────────┐──────────────┐
    │ chrA:101-200 │ chrA:201-300 │ chrA:301-400 │ chrA:401-500 │ (Target genome)
    └──────────────┴──────────────┴──────────────┘──────────────┘
           +               -             -              +   (Alignment direction)
    ┌──────────────┬──────────────┬──────────────┐──────────────┐
    │ chrB:101-200 │ chrB:201-300 │ chrB:301-400 │ chrB:401-500 │  (Query genome)
    └──────────────┴──────────────┴──────────────┘──────────────┘

## See also

Other Flagging functions:
[`flagAll()`](https://oist.github.io/GenomicBreaks/reference/flagAll.md),
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagLongShort()`](https://oist.github.io/GenomicBreaks/reference/flagLongShort.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md)

Other Inversion functions:
[`filterDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/filterDoubleInversions.md),
[`filterInversions()`](https://oist.github.io/GenomicBreaks/reference/filterInversions.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`leftInversionGaps()`](https://oist.github.io/GenomicBreaks/reference/leftInversionGaps.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
[`showInversions()`](https://oist.github.io/GenomicBreaks/reference/showInversions.md)

Other Structural variants:
[`StructuralVariants`](https://oist.github.io/GenomicBreaks/reference/StructuralVariants.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
[`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)

## Examples

``` r
flagTwinInversions(exampleTwinInversions)
#> GBreaks object with 4 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query       twi
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <logical>
#>   [1]     chrA   100-190      + | chrB:100-190      TRUE
#>   [2]     chrA   200-290      - | chrB:200-290     FALSE
#>   [3]     chrA   300-390      - | chrB:300-390     FALSE
#>   [4]     chrA   400-490      + | chrB:400-490     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
plotApairOfChrs(exampleTwinInversions)


flagTwinInversions(exampleTwinInversions |> reverse() |> sort(ignore.strand = TRUE))
#> GBreaks object with 4 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query       twi
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <logical>
#>   [1]     chrA   111-201      - | chrB:400-490      TRUE
#>   [2]     chrA   211-301      + | chrB:300-390     FALSE
#>   [3]     chrA   311-401      + | chrB:200-290     FALSE
#>   [4]     chrA   411-501      - | chrB:100-190     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
plotApairOfChrs(exampleTwinInversions |> reverse())

```
