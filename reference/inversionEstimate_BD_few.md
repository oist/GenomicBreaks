# Inversion Estimate - option 'few'

It finds the number of inversions whose expected number of cycles is
closest to the observed number of cycles.

## Usage

``` r
inversionEstimate_BD_few(n, obs_nb_cycles)
```

## Arguments

- n:

  An integer: The number of markers. A marker represents an aligned
  block, a gene, or any other genomic region being rearranged.

- obs_nb_cycles:

  An integer: The observed number of cycles in the breakpoint graph of
  two genomes.

## Details

This function is an alternative to the function
"inversionEstimate_BD_many", which uses the function
"inversionEstimate_BD_setup".

This function should be preferred when there are few inversion estimates
to be computed.

## Author

Priscila Biller
