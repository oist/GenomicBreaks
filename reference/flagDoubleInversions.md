# Flag double inversions

Flag ranges that start a quintuplet in which two overlapping inversions
moved the three central ranges, and the two flanking ranges did not
move.

## Usage

``` r
flagDoubleInversions(gb, details = FALSE)
```

## Arguments

- gb:

  A
  [`GBreaks()`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- details:

  Returns more metadata columns if `TRUE`.

## Value

Returns the `GBreaks` object with an extra `Dbl` metadata column.

## Details

With the five *target* ranges named `X A B C Y` in their order on the
*target* genome, and lower case meaning the minus strand, the *query*
genome reads one of the following:

    X b c A Y    ABC -> baC -> bcA
    X C a b Y    ABC -> Acb -> Cab
    y B A c x    cba -> BCa -> BAc
    y a C B x    cba -> cAB -> aCB

The last two are the first two read on the other strand of the *query*
genome. The flanks `X` and `Y` must be on the strand shown: with flanks
on the opposite strand, the same central triplet requires three
inversions. The flag is set on range `X`, as
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md)
sets its flag on the left flank of an inversion. The five ranges must be
on the same *target* sequence, and their strand must be known (not `*`).

## See also

Other Flagging functions:
[`flagAll()`](https://oist.github.io/GenomicBreaks/reference/flagAll.md),
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagLongShort()`](https://oist.github.io/GenomicBreaks/reference/flagLongShort.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
[`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md)

Other Inversion functions:
[`filterDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/filterDoubleInversions.md),
[`filterInversions()`](https://oist.github.io/GenomicBreaks/reference/filterInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`leftInversionGaps()`](https://oist.github.io/GenomicBreaks/reference/leftInversionGaps.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
[`showInversions()`](https://oist.github.io/GenomicBreaks/reference/showInversions.md)

Other Structural variants:
[`StructuralVariants`](https://oist.github.io/GenomicBreaks/reference/StructuralVariants.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
[`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md),
[`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)

## Author

Charles Plessy

## Examples

``` r
# Start colinear.  Lower case meands minus strand
z0 <- GBreaks(target = GRanges(c(A="T:10-15:+", B="T:20-25:+", C="T:30-35:+")),
              query  = GRanges(c(A="Q:10-15",   B="Q:20-25",   C="Q:30-35")))

# Swap coordinates of A and B on the query and flip strands  ABC -> baC
z1 <- GBreaks(target = GRanges(c(A="T:10-15:-", B="T:20-25:-", C="T:30-35:+")),
              query  = GRanges(c(a="Q:20-25",   b="Q:10-15",   C="Q:30-35")))

# Now query order is b - a - C.  Swap a and C and flip strands baC -> bcA
z2 <- GBreaks(target = GRanges(c(A="T:10-15:+", B="T:20-25:-", C="T:30-35:-")),
              query  = GRanges(c(A="Q:30-35",   b="Q:10-15",   c="Q:20-25")))

# Altogether, there are:
# ABC -> baC -> bcA
# ABC -> Acb -> Cab
# cba -> BCa -> BAc
# cba -> cAB -> aCB


# z2 has same topoloty as package example
plotApairOfChrs(z2)

plotApairOfChrs(exampleDoubleInversion1)


flagDoubleInversions(exampleDoubleInversion1 )
#> GBreaks object with 5 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query       Dbl
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <logical>
#>   [1]     chrA   100-190      + | chrB:100-190      TRUE
#>   [2]     chrA   200-290      + | chrB:400-490     FALSE
#>   [3]     chrA   300-390      - | chrB:200-290     FALSE
#>   [4]     chrA   400-490      - | chrB:300-390     FALSE
#>   [5]     chrA   500-590      + | chrB:500-590     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome

exampleDoubleInversion1Rev <- reverse(exampleDoubleInversion1) |> sort(ignore.strand = TRUE)
exampleDoubleInversion1Rev[2:4] |> plotApairOfChrs()

```
