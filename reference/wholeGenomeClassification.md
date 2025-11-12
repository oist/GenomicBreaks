# Whole-genome classification object

Classifies aligned genomic regions as *isolated* or *collinear*, and
unaligned regions as *breakpoint* or *bridge* regions. The *bridge*
regions are always flanked by *collinear alignments* and the *isolated
alignments* are always flanked by *breakpoint regions*.

## Usage

``` r
wholeGenomeClassification(gb, coa = coalesce_contigs(gb), ends = FALSE)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object representing a one-to-one whole genome alignment.

- coa:

  The coalesced one-to-one alignment. If not provided, it will be\\
  computed on-the-fly with the
  [`coalesce_contigs`](https://oist.github.io/GenomicBreaks/reference/coalesce_contigs.md)
  function.

- ends:

  Add a a *end region* type for the extremities of the *sequence
  features* not covered by the original alignment.

## Value

A
[`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object representing the *target* genome, inheriting its *sequence
information*
([`GenomeInfoDb::Seqinfo`](https://rdrr.io/pkg/Seqinfo/man/Seqinfo-class.html)).
The class of each region is indicated by a factor in the `type` metadata
column.

## Details

The strand of *bridge regions* is the one of their flanking *collinear
alignment* regions. *Breakpoint* and *end* regions are unstranded.

## See also

Other Reducing functions:
[`cleanGaps()`](https://oist.github.io/GenomicBreaks/reference/cleanGaps.md),
[`get_bps()`](https://oist.github.io/GenomicBreaks/reference/get_bps.md)

## Author

Charles Plessy

Michael Mansfield

## Examples

``` r
exampleColinear5 |> wholeGenomeClassification(ends = TRUE)
#> GRanges object with 11 ranges and 1 metadata column:
#>        seqnames    ranges strand |                type
#>           <Rle> <IRanges>  <Rle> |            <factor>
#>    [1]     chrA      1-99      * | end region         
#>    [2]     chrA   100-190      + | collinear alignment
#>    [3]     chrA   191-199      + | bridge region      
#>    [4]     chrA   200-290      + | collinear alignment
#>    [5]     chrA   291-299      + | bridge region      
#>    [6]     chrA   300-390      + | collinear alignment
#>    [7]     chrA   391-399      + | bridge region      
#>    [8]     chrA   400-490      + | collinear alignment
#>    [9]     chrA   491-499      + | bridge region      
#>   [10]     chrA   500-590      + | collinear alignment
#>   [11]     chrA   591-600      * | end region         
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
