# Inversion Estimate - setup for option 'many'

The setup step computes the expected number of cycles after k edges for
various values of k.

## Usage

``` r
inversionEstimate_BD_setup(n)
```

## Arguments

- n:

  An integer: The number of markers. A marker represents an aligned
  block, a gene, or any other genomic region being rearranged.

## Details

The expected number of cycles can take any value between 1 and n (the
number of markers). For each of these values, the method computes at
least one k.

Notice that this is the most time-consuming step in the estimate. Since
it depends only on the number n of markers, it is recommended that it be
computed once for all permutations of n markers, saving time.

## Author

Priscila Biller
