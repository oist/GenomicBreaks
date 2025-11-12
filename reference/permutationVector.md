# Permutation vector

Represent *GenomicBreaks* objects as a *permutation vector*.

## Usage

``` r
permutationVector(gb)
```

## Arguments

- gb:

  A
  [`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
  object.

## Value

A vector of non-zero integers representing the position of the genomic
content in the `GBreaks` object.

## Details

In a `GBreaks` object properly sorted by its *target* ranges, the
*query* ranges …
