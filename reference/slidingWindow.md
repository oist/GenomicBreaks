# Sliding window

This function is used to subset a `GBreaks` object into a list of
`GBreaks` using a sliding window.

## Usage

``` r
slidingWindow(
  gb,
  factor = 100,
  windowSize = NULL,
  stepSize = NULL,
  merged = FALSE,
  cut = TRUE
)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

- factor:

  Number of windows.

- windowSize:

  The size of the window.

- stepSize:

  The size of the step.

- merged:

  set TRUE to output a `GBreaks` object

- cut:

  whether cut the alignments out of the window

## Value

Returns a
[`GenomicRanges::GRangesList`](https://rdrr.io/pkg/GenomicRanges/man/GRangesList-class.html)
containing one `GBreaks` object for each window, or a single `GBreaks`
object if `merged` is `TRUE`.

## Details

When `cut` is `TRUE`, an alignment that crosses a window boundary is cut
at that boundary on the *target* genome. Its *query* range is cut by the
same proportion, on the side given by the strand of the *target* range,
and its score, if any, is scaled by that proportion.

## Unstranded ranges (experimental)

The *query* ranges of alignments on the `*` strand are not cut.

## See also

Other Similarity indexes:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`JC69_distance_allseq()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance_allseq.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`K80_gap_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_gap_distance.md),
[`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md),
[`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md),
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`breakpointGraphProperties()`](https://oist.github.io/GenomicBreaks/reference/breakpointGraphProperties.md),
[`correlation_index()`](https://oist.github.io/GenomicBreaks/reference/correlation_index.md),
[`inversionDistance()`](https://oist.github.io/GenomicBreaks/reference/inversionDistance.md),
[`inversionEstimate_BD()`](https://oist.github.io/GenomicBreaks/reference/inversionEstimate_BD.md),
[`karyotype_index()`](https://oist.github.io/GenomicBreaks/reference/karyotype_index.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md),
[`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md),
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md),
[`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)

## Author

Zhang Kun

Claude Opus 5.5 (Anthropic), prompted by Charles Plessy

## Examples

``` r
exampleInversion
#> GBreaks object with 3 ranges and 1 metadata column:
#>       seqnames    ranges strand |        query
#>          <Rle> <IRanges>  <Rle> |    <GRanges>
#>   [1]     chrA   100-190      + | chrB:100-190
#>   [2]     chrA   210-291      - | chrB:210-291
#>   [3]     chrA   301-400      + | chrB:301-400
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
slided_gb <- slidingWindow(exampleInversion, windowSize = 250, stepSize = 125)
slided_gb
#> GRangesList object of length 3:
#> $`1`
#> GBreaks object with 2 ranges and 5 metadata columns:
#>       seqnames    ranges strand |        query window_id window_seqname
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <integer>    <character>
#>   [1]     chrA   100-190      + | chrB:100-190         1           chrA
#>   [2]     chrA   210-250      - | chrB:251-291         1           chrA
#>       window_start window_end
#>          <integer>  <integer>
#>   [1]            1        250
#>   [2]            1        250
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
#> 
#> $`2`
#> GBreaks object with 3 ranges and 5 metadata columns:
#>       seqnames    ranges strand |        query window_id window_seqname
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <integer>    <character>
#>   [1]     chrA   126-190      + | chrB:126-190         2           chrA
#>   [2]     chrA   210-291      - | chrB:210-291         2           chrA
#>   [3]     chrA   301-375      + | chrB:301-375         2           chrA
#>       window_start window_end
#>          <integer>  <integer>
#>   [1]          126        375
#>   [2]          126        375
#>   [3]          126        375
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
#> 
#> $`3`
#> GBreaks object with 2 ranges and 5 metadata columns:
#>       seqnames    ranges strand |        query window_id window_seqname
#>          <Rle> <IRanges>  <Rle> |    <GRanges> <integer>    <character>
#>   [1]     chrA   251-291      - | chrB:210-250         3           chrA
#>   [2]     chrA   301-400      + | chrB:301-400         3           chrA
#>       window_start window_end
#>          <integer>  <integer>
#>   [1]          251        500
#>   [2]          251        500
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
#> 
sapply(slided_gb, strand_randomisation_index)
#>         1         2         3 
#> 0.3787879 0.2612613 0.4184397 
```
