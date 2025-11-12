# Flag inversions

Flag ranges that start a triplet that would be colinear if the central
pair were aligned to the opposite strand.

## Usage

``` r
flagInversions(gb, tol = Inf)
```

## Arguments

- gb:

  A
  [`GBreaks()`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- tol:

  Tolerance window for the distance between two ranges.

## Value

Returns the `GBreaks` object with an extra `inv` metadata column.

## Details

Here is a trivial example of an inversion.

    ┌──────────────┬──────────────┬──────────────┐
    │ chrA:101-200 │ chrA:201-300 │ chrA:301-400 │ (Target genome)
    └──────────────┴──────────────┴──────────────┘
           +               -             +         (Alignment direction)
    ┌──────────────┬──────────────┬──────────────┐
    │ chrB:101-200 │ chrB:201-300 │ chrB:301-400 │ (Query genome)
    └──────────────┴──────────────┴──────────────┘

## See also

Other Flagging functions:
[`flagAll()`](https://oist.github.io/GenomicBreaks/reference/flagAll.md),
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagLongShort()`](https://oist.github.io/GenomicBreaks/reference/flagLongShort.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
[`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md)

Other Inversion functions:
[`filterDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/filterDoubleInversions.md),
[`filterInversions()`](https://oist.github.io/GenomicBreaks/reference/filterInversions.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`leftInversionGaps()`](https://oist.github.io/GenomicBreaks/reference/leftInversionGaps.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
[`showInversions()`](https://oist.github.io/GenomicBreaks/reference/showInversions.md)

Other Structural variants:
[`StructuralVariants`](https://oist.github.io/GenomicBreaks/reference/StructuralVariants.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
[`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md),
[`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)

## Examples

``` r
flagInversions(exampleInversion)
#> GBreaks object with 3 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query   inv
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle>
#>   [1]     chrA   100-190      + | chrB:100-190  TRUE
#>   [2]     chrA   210-291      - | chrB:210-291 FALSE
#>   [3]     chrA   301-400      + | chrB:301-400 FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
flagInversions(exampleInversion, tol = 19)
#> GBreaks object with 3 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query   inv
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle>
#>   [1]     chrA   100-190      + | chrB:100-190 FALSE
#>   [2]     chrA   210-291      - | chrB:210-291 FALSE
#>   [3]     chrA   301-400      + | chrB:301-400 FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
plotApairOfChrs(exampleInversion)


flagInversions(exampleInversion |> reverse() |> sort(ignore.strand = TRUE))
#> GBreaks object with 3 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query   inv
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle>
#>   [1]     chrA   201-300      - | chrB:301-400  TRUE
#>   [2]     chrA   310-391      + | chrB:210-291 FALSE
#>   [3]     chrA   411-501      - | chrB:100-190 FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
plotApairOfChrs(exampleInversion |> reverse())

```
