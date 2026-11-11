# Flag translocations

Flag the two places where a range that moved between the two genomes is
seen in the order of the *target* genome: the *payload region*, where
the moved range is found on the *target* genome, and the *site*, where
it is found on the *query* genome.

## Usage

``` r
flagTranslocations(gb, tol = Inf)
```

## Arguments

- gb:

  A
  [`GBreaks()`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- tol:

  Ranges separated by more than `tol` unaligned bases are not considered
  adjacent: `U`, `A` and `V` on the *target* genome and `U` and `V` on
  the *query* genome for payload regions; `X` and `Y` on the *target*
  genome and `X`, `A` and `Y` on the *query* genome for sites.

## Value

Returns the `GBreaks` object with two extra metadata columns, `tra` and
`trs`, flagging the payload regions and the sites.

## Details

A range `A` moved between the two genomes when, for instance:

    Target genome:  X Y . . . U A V
    Query genome:   X A Y . . . U V

- The *payload region* is a triplet of consecutive ranges `U A V` on one
  *target* sequence, such that `U` and `V` are on the same strand and
  would be colinear if `A` were removed. The `tra` flag is set on `U`.

- The *site* is a pair of consecutive ranges `X Y` on one *target*
  sequence, on the same strand, such that `A` is found between them on
  the *query* genome, and that `X A Y` would be colinear if `A` were on
  the same strand as `X` and `Y`. The `trs` flag is set on `X`.

`A` can be on either strand, and on any *target* sequence. When one of
the two genomes has no context for `A`, for instance when it is alone on
its sequence, only one of the two flags is set. When both are set, they
designate the same range `A`. The function therefore detects
transpositions of single ranges, whether or not they change sequence.

Here is a trivial example of a translocation, where the *query* genome
has no context for the moved range, so that only the payload region is
flagged:

    ┌──────────────┬──────────────┬──────────────┐
    │ chrA:101-200 │ chrA:201-300 │ chrA:301-400 │ (Target genome)
    └──────────────┴──────────────┴──────────────┘
           +               +             +         (Alignment direction)
    ┌──────────────┬──────────────┬──────────────┐
    │ chrB:101-200 │ chrC:201-300 │ chrB:301-400 │ (Query genome)
    └──────────────┴──────────────┴──────────────┘

## Note

Before version 0.25.0, the site was flagged only with the `both` option,
by a `tra` flag on the *query* ranges, set on the first range of the
site in *query* order.

## Unstranded ranges (experimental)

Ranges on the `*` strand are considered adjacent in either direction.
When orientation is unknown, moving one range over a segment can not be
told apart from moving the segment over the range, so that more ranges
are flagged as moved: for instance, with the *target* genome `X Y U A V`
and the *query* genome `X A Y U V`, `Y` and `U` are flagged in addition
to `A`. This behaviour may change in future versions.

## See also

Other Flagging functions:
[`flagAll()`](https://oist.github.io/GenomicBreaks/reference/flagAll.md),
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagLongShort()`](https://oist.github.io/GenomicBreaks/reference/flagLongShort.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md)

Other Translocation functions:
[`filterTranslocations()`](https://oist.github.io/GenomicBreaks/reference/filterTranslocations.md),
[`removeTranslocations()`](https://oist.github.io/GenomicBreaks/reference/removeTranslocations.md),
[`showTranslocations()`](https://oist.github.io/GenomicBreaks/reference/showTranslocations.md)

Other Structural variants:
[`StructuralVariants`](https://oist.github.io/GenomicBreaks/reference/StructuralVariants.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagPairs()`](https://oist.github.io/GenomicBreaks/reference/flagPairs.md),
[`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md),
[`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)

## Examples

``` r
flagTranslocations(exampleTranslocation)
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query   tra   trs
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   [1]     chrA   100-200      + | chrB:100-200  TRUE FALSE
#>   [2]     chrA   201-300      + | chrC:201-300 FALSE FALSE
#>   [3]     chrA   301-400      + | chrB:301-400 FALSE FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
# The ranges U and V are separated by 100 bases on the query genome.
flagTranslocations(exampleTranslocation, tol = 100)
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query   tra   trs
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   [1]     chrA   100-200      + | chrB:100-200  TRUE FALSE
#>   [2]     chrA   201-300      + | chrC:201-300 FALSE FALSE
#>   [3]     chrA   301-400      + | chrB:301-400 FALSE FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
flagTranslocations(exampleTranslocation, tol = 99)
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query   tra   trs
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   [1]     chrA   100-200      + | chrB:100-200 FALSE FALSE
#>   [2]     chrA   201-300      + | chrC:201-300 FALSE FALSE
#>   [3]     chrA   301-400      + | chrB:301-400 FALSE FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
flagTranslocations(exampleTranslocation2)
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query   tra   trs
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   [1]     chrA   100-200      + | chrB:100-200  TRUE FALSE
#>   [2]     chrA   201-300      - | chrC:201-300 FALSE FALSE
#>   [3]     chrA   301-400      + | chrB:301-400 FALSE FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
plotApairOfChrs(exampleTranslocation)


# When the genomes are swapped, the moved range has no context on the
# target genome, and only the site is flagged.
flagTranslocations(swap(exampleTranslocation, sort = TRUE))
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query   tra   trs
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   [1]     chrB   100-200      + | chrA:100-200 FALSE  TRUE
#>   [2]     chrB   301-400      + | chrA:301-400 FALSE FALSE
#>   [3]     chrC   201-300      + | chrA:201-300 FALSE FALSE
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome

# A move along one sequence is seen at both places.
gb <- permutationVectorToGBreaks(c(1, 3, 4, 2, 5))
flagTranslocations(gb)
#> GBreaks object with 5 ranges and 3 metadata columns:
#>       seqnames    ranges strand |         query   tra   trs
#>          <Rle> <IRanges>  <Rle> |     <GRanges> <Rle> <Rle>
#>   [1]    chrT1   100-189      + | chrQ1:100-189 FALSE  TRUE
#>   [2]    chrT1   200-289      + | chrQ1:300-389 FALSE FALSE
#>   [3]    chrT1   300-389      + | chrQ1:400-489  TRUE FALSE
#>   [4]    chrT1   400-489      + | chrQ1:200-289 FALSE FALSE
#>   [5]    chrT1   500-589      + | chrQ1:500-589 FALSE FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome

flagTranslocations(exampleDeletion)
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query   tra   trs
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   [1]     chrA   100-200      + | chrB:100-200  TRUE FALSE
#>   [2]     chrA   201-300      + | chrC:401-500 FALSE FALSE
#>   [3]     chrA   301-400      + | chrB:201-300 FALSE FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
flagTranslocations(exampleInsertion)
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query   tra   trs
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   [1]     chrA   100-200      + | chrB:100-200 FALSE  TRUE
#>   [2]     chrA   201-300      + | chrB:301-400 FALSE FALSE
#>   [3]     chrC   401-500      + | chrB:201-300 FALSE FALSE
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
flagTranslocations(exampleInsertion |> swap(sort = TRUE))
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query   tra   trs
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   [1]     chrB   100-200      + | chrA:100-200  TRUE FALSE
#>   [2]     chrB   201-300      + | chrC:401-500 FALSE FALSE
#>   [3]     chrB   301-400      + | chrA:201-300 FALSE FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
flagTranslocations(sort(reverse(exampleDeletion), ignore.strand = TRUE))
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query   tra   trs
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   [1]     chrA   201-300      - | chrB:201-300  TRUE FALSE
#>   [2]     chrA   301-400      - | chrC:401-500 FALSE FALSE
#>   [3]     chrA   401-501      - | chrB:100-200 FALSE FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
flagTranslocations(exampleInversion)
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query   tra   trs
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   [1]     chrA   100-190      + | chrB:100-190 FALSE FALSE
#>   [2]     chrA   210-291      - | chrB:210-291 FALSE FALSE
#>   [3]     chrA   301-400      + | chrB:301-400 FALSE FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
flagTranslocations(exampleColinear3)
#> GBreaks object with 3 ranges and 3 metadata columns:
#>       seqnames    ranges strand |        query   tra   trs
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <Rle> <Rle>
#>   [1]     chrA   100-200      + | chrB:100-200 FALSE FALSE
#>   [2]     chrA   201-300      + | chrB:201-300 FALSE FALSE
#>   [3]     chrA   301-400      + | chrB:301-400 FALSE FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
