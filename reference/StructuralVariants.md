# Structural Variants

Concept examples of structural variants for test and demonstration
purpose

## Format

[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
objects described in more details in the *Structural Variants* vignette,
available with the command
[`vignette("StructuralVariants", package = "GenomicBreaks")`](https://oist.github.io/GenomicBreaks/articles/StructuralVariants.md).

## Details

The variants represented here are the ones that can be represented in a
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
object. Short ones that are entirely contained in an aligned region, for
instance as a gap, are not covered

## See also

Other Structural variants:
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
[`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)

## Author

Charles Plessy

## Examples

``` r
exampleColinear
#> GBreaks object with 2 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-150      + | chrB:100-150
#>   [2]     chrA   251-300      + | chrB:251-300
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
exampleNotColinear
#> GBreaks object with 2 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-150      + | chrB:201-251
#>   [2]     chrA   251-300      + |  chrB:50-100
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
exampleColinear3
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-200      + | chrB:100-200
#>   [2]     chrA   201-300      + | chrB:201-300
#>   [3]     chrA   301-400      + | chrB:301-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
exampleColinear5
#> GBreaks object with 5 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-190      + | chrB:100-190
#>   [2]     chrA   200-290      + | chrB:200-290
#>   [3]     chrA   300-390      + | chrB:300-390
#>   [4]     chrA   400-490      + | chrB:400-490
#>   [5]     chrA   500-590      + | chrB:500-590
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
exampleInversion
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-190      + | chrB:100-190
#>   [2]     chrA   210-291      - | chrB:210-291
#>   [3]     chrA   301-400      + | chrB:301-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
exampleInversion5uncollapsed
#> GBreaks object with 5 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-190      + | chrB:100-190
#>   [2]     chrA   200-290      - | chrB:300-390
#>   [3]     chrA   300-390      - | chrB:200-290
#>   [4]     chrA   400-490      + | chrB:400-490
#>   [5]     chrA   500-590      + | chrB:500-590
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
exampleNestedInversions
#> GBreaks object with 5 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-190      + | chrA:100-190
#>   [2]     chrA   200-290      - | chrA:400-490
#>   [3]     chrA   300-390      + | chrA:300-390
#>   [4]     chrA   400-490      - | chrA:200-290
#>   [5]     chrA   500-590      + | chrA:500-590
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
exampleDoubleInversion1
#> GBreaks object with 5 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-190      + | chrB:100-190
#>   [2]     chrA   200-290      + | chrB:400-490
#>   [3]     chrA   300-390      - | chrB:200-290
#>   [4]     chrA   400-490      - | chrB:300-390
#>   [5]     chrA   500-590      + | chrB:500-590
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
exampleTranslocation
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-200      + | chrB:100-200
#>   [2]     chrA   201-300      + | chrC:201-300
#>   [3]     chrA   301-400      + | chrB:301-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
exampleTranslocation2
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-200      + | chrB:100-200
#>   [2]     chrA   201-300      - | chrC:201-300
#>   [3]     chrA   301-400      + | chrB:301-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
