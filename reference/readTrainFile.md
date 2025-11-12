# Extract values from `last-train` parameter file

This function is used to read the train file from LAST and extract the
final parameters.

## Usage

``` r
readTrainFile(input_file)
```

## Arguments

- input_file:

  A string of the path to the train file from LAST.

## Value

Returns a `list` of parameters from LAST train file.
"PercentSimilarity", "PercentSimilarityNogaps", "mean_delete_size",
"mean_insert_size", "substitution_percent_identity", "probability",
"matchProb", "delOpenProb", "insOpenProb", "delExtendProb",
"insExtendProb", "endProb"

## See also

Other Data loading functions:
[`load_genomic_breaks()`](https://oist.github.io/GenomicBreaks/reference/load_genomic_breaks.md),
[`readMAF()`](https://oist.github.io/GenomicBreaks/reference/readMAF.md)

## Author

Zikun Yang

## Examples

``` r
readTrainFile(system.file("extdata/example.train", package = "GenomicBreaks"))
#> $mean_delete_size
#> [1] 6.25287
#> 
#> $mean_insert_size
#> [1] 6.80408
#> 
#> $substitution_percent_identity
#> [1] 98.4025
#> 
#> $probability_A_A
#> [1] 0.289718
#> 
#> $probability_A_C
#> [1] 0.000727394
#> 
#> $probability_A_G
#> [1] 0.00253215
#> 
#> $probability_A_T
#> [1] 0.000611873
#> 
#> $probability_C_A
#> [1] 0.000689766
#> 
#> $probability_C_C
#> [1] 0.202294
#> 
#> $probability_C_G
#> [1] 0.000810395
#> 
#> $probability_C_T
#> [1] 0.00261617
#> 
#> $probability_G_A
#> [1] 0.00261617
#> 
#> $probability_G_C
#> [1] 0.000810395
#> 
#> $probability_G_G
#> [1] 0.202294
#> 
#> $probability_G_T
#> [1] 0.000689766
#> 
#> $probability_T_A
#> [1] 0.000611873
#> 
#> $probability_T_C
#> [1] 0.00253215
#> 
#> $probability_T_G
#> [1] 0.000727394
#> 
#> $probability_T_T
#> [1] 0.289718
#> 
#> $matchProb
#> [1] 0.997811
#> 
#> $delOpenProb
#> [1] 0.000895537
#> 
#> $insOpenProb
#> [1] 0.000766701
#> 
#> $delExtendProb
#> [1] 0.795939
#> 
#> $insExtendProb
#> [1] 0.795939
#> 
#> $endProb
#> [1] 0.999737
#> 
#> $probability_matrix
#>             A           C           G           T
#> A 0.289718000 0.000727394 0.002532150 0.000611873
#> C 0.000689766 0.202294000 0.000810395 0.002616170
#> G 0.002616170 0.000810395 0.202294000 0.000689766
#> T 0.000611873 0.002532150 0.000727394 0.289718000
#> 
```
