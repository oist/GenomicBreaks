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

Returns a list of `GBreaks` objects, one for each window.

## See also

Other Similarity indexes:
[`F81_distance()`](https://oist.github.io/GenomicBreaks/reference/F81_distance.md),
[`GOC()`](https://oist.github.io/GenomicBreaks/reference/GOC.md),
[`HKY85_distance()`](https://oist.github.io/GenomicBreaks/reference/HKY85_distance.md),
[`JC69_distance()`](https://oist.github.io/GenomicBreaks/reference/JC69_distance.md),
[`K80_distance()`](https://oist.github.io/GenomicBreaks/reference/K80_distance.md),
[`P_distance()`](https://oist.github.io/GenomicBreaks/reference/P_distance.md),
[`T92_distance()`](https://oist.github.io/GenomicBreaks/reference/T92_distance.md),
[`TN93_distance()`](https://oist.github.io/GenomicBreaks/reference/TN93_distance.md),
[`correlation_index()`](https://oist.github.io/GenomicBreaks/reference/correlation_index.md),
[`karyotype_index()`](https://oist.github.io/GenomicBreaks/reference/karyotype_index.md),
[`logDet_distance()`](https://oist.github.io/GenomicBreaks/reference/logDet_distance.md),
[`strand_randomisation_index()`](https://oist.github.io/GenomicBreaks/reference/strand_randomisation_index.md),
[`synteny_index()`](https://oist.github.io/GenomicBreaks/reference/synteny_index.md),
[`tau_index()`](https://oist.github.io/GenomicBreaks/reference/tau_index.md)

## Author

Zhang Kun

## Examples

``` r
exdata_Sac <- system.file("extdata/SacCer3__SacPar.gff3.gz", package = "GenomicBreaks")
gb <- load_genomic_breaks(exdata_Sac, BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae)
slided_gb <- slidingWindow(gb, windowSize = 1e5, stepSize = 5e4)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
slided_gb[[1]]
#> GBreaks object with 6 ranges and 6 metadata columns:
#>       seqnames       ranges strand |     score                    query
#>          <Rle>    <IRanges>  <Rle> | <numeric>                <GRanges>
#>   [1]     chrI   5860-10010      + |     10609   NC_047487.1:8723-12716
#>   [2]     chrI  11157-11782      + |      1253  NC_047487.1:13035-13680
#>   [3]     chrI  25371-26528      + |      3651  NC_047487.1:15380-16537
#>   [4]     chrI  26849-29699      + |      6130  NC_047487.1:17263-20175
#>   [5]     chrI  29937-30607      + |       767  NC_047487.1:20203-20866
#>   [6]     chrI 31155-100000      + |    298623 NC_047487.1:21379-128635
#>       window_id window_seqname window_start window_end
#>       <integer>    <character>    <integer>  <integer>
#>   [1]         1           chrI            1     100000
#>   [2]         1           chrI            1     100000
#>   [3]         1           chrI            1     100000
#>   [4]         1           chrI            1     100000
#>   [5]         1           chrI            1     100000
#>   [6]         1           chrI            1     100000
#>   -------
#>   seqinfo: 17 sequences (1 circular) from sacCer3 genome
```
