# Divergence time - TimeTree

Retrieves the divergence time of two species via the TimeTree API.

## Usage

``` r
getDivergenceTime(taxid1, taxid2)
```

## Arguments

- taxid1:

  An integer: The NCBI id of one taxon

- taxid2:

  An integer: The NCBI id of the other taxon

## Value

A list in which the first item is the median divergence time
(`divtimeMedian`), and the last two items are the divergence time
confidence interval (`divtimeCI_low` and `divtimeCI_high`). If the
divergence time is not found in the TimeTree database, it returns `NA`.

## Author

Priscila Biller

## Examples

``` r
if (FALSE) { # \dontrun{
getDivergenceTime(9555,9601)} # }
```
