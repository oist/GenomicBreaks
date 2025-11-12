# Flag successive entries of a `GBreaks` object

Scans a sorted
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
object and issues a flag describing the relation between the current
entry and the next entry.

## Usage

``` r
flagPairs(gb)
```

## Arguments

- gb:

  A
  [`GBreaks()`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

## Value

Returns the `GBreaks` object with an extra `pairs` metadata column. This
`pairs` column is a factor of all flags, so that `table(gb$pairs)`
always returns a vector of the same length, reporting `0` when
necessary. This allows to easily aggregate results from multiple
objects.

## Details

Because the object is sorted, the next entry on the *target* genome is
by definition following the current one unless the end of the sequence
feature (contig, scaffold, …) is reached. The two *target* ranges can be
on the same strand or on opposite strands. No difference is made between
`+/-` and `-/+` orientations because it will be arbitrary unless the
sequence features of both genomes are oriented, as in the case of
comparing chromosomal assemblies of related species.

Depending on whether the *query* ranges precede or follow each other, or
are not next to each other, and depending on the strand of the
alignments, a flag is issued on the current entry. `Col` signals
colinearity with the next entry, `Inv` an inversion on either entry,
`Flp` signals that the order of the entries is as if one hand jumped
over the other one. When the *query* ranges are not next to each other,
a `Scr` flag, for *scrambled*, is issued. Lastly, `Bnd` (for *boundary*)
signals that there is no pair to analyse because the current entry is
the last one for the current sequence level on the *target* genome.

## See also

Other Flagging functions:
[`flagAll()`](https://oist.github.io/GenomicBreaks/reference/flagAll.md),
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagLongShort()`](https://oist.github.io/GenomicBreaks/reference/flagLongShort.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
[`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md)

Other Inversion functions:
[`filterDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/filterDoubleInversions.md),
[`filterInversions()`](https://oist.github.io/GenomicBreaks/reference/filterInversions.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md),
[`flipInversions()`](https://oist.github.io/GenomicBreaks/reference/flipInversions.md),
[`leftInversionGaps()`](https://oist.github.io/GenomicBreaks/reference/leftInversionGaps.md),
[`removeDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/removeDoubleInversions.md),
[`removeInversions()`](https://oist.github.io/GenomicBreaks/reference/removeInversions.md),
[`showInversions()`](https://oist.github.io/GenomicBreaks/reference/showInversions.md)

Other Colinearity functions:
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`bridgeRegions()`](https://oist.github.io/GenomicBreaks/reference/bridgeRegions.md),
[`chain_contigs()`](https://oist.github.io/GenomicBreaks/reference/chain_contigs.md),
[`coalesce_contigs()`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md),
[`dist2next()`](https://oist.github.io/GenomicBreaks/reference/dist2next.md),
[`filterColinearRegions()`](https://oist.github.io/GenomicBreaks/reference/filterColinearRegions.md),
[`flagColinearAlignments()`](https://oist.github.io/GenomicBreaks/reference/flagColinearAlignments.md)

Other Structural variants:
[`StructuralVariants`](https://oist.github.io/GenomicBreaks/reference/StructuralVariants.md),
[`flagDoubleInversions()`](https://oist.github.io/GenomicBreaks/reference/flagDoubleInversions.md),
[`flagInversions()`](https://oist.github.io/GenomicBreaks/reference/flagInversions.md),
[`flagTranslocations()`](https://oist.github.io/GenomicBreaks/reference/flagTranslocations.md),
[`flagTwinInversions()`](https://oist.github.io/GenomicBreaks/reference/flagTwinInversions.md),
[`plotApairOfChrs()`](https://oist.github.io/GenomicBreaks/reference/plotApairOfChrs.md)

## Examples

``` r
flagPairs(exampleInversion)
#> GBreaks object with 3 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query    pairs
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <factor>
#>   [1]     chrA   100-190      + | chrB:100-190      Inv
#>   [2]     chrA   210-291      - | chrB:210-291      Inv
#>   [3]     chrA   301-400      + | chrB:301-400      Bnd
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
flagPairs(exampleTranslocation)
#> GBreaks object with 3 ranges and 2 metadata columns:
#>       seqnames    ranges strand |        query    pairs
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <factor>
#>   [1]     chrA   100-200      + | chrB:100-200      Scr
#>   [2]     chrA   201-300      + | chrC:201-300      Scr
#>   [3]     chrA   301-400      + | chrB:301-400      Bnd
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome

# How the pair classes are defined:
GenomicBreaks:::allPossiblePairClasses
#>    Qnext Strnd Snext pair     paste
#> 1   next     +     +  Col  next + +
#> 2   prev     +     +  Flp  prev + +
#> 3  other     +     +  Scr other + +
#> 4   next     -     +  Inv  next - +
#> 5   prev     -     +  Inv  prev - +
#> 6  other     -     +  Scr other - +
#> 7   next     +     -  Inv  next + -
#> 8   prev     +     -  Inv  prev + -
#> 9  other     +     -  Scr other + -
#> 10  next     -     -  Flp  next - -
#> 11  prev     -     -  Col  prev - -
#> 12 other     -     -  Scr other - -
```
