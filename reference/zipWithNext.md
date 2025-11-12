# Zip elements with the next ones

Zip an object with itself after inserting a variable number of steps.
For instance, with a step of 1, *a* *b* *c* would become *a,b*, *b,c*,
*c,NA*.

## Usage

``` r
zipWithNext(x, step = 1)
```

## Arguments

- x:

  Object with a vector structure.

- step:

  Remove the first `step` entries from the object before zipping it with
  itself.

## Value

A
[`S4Vectors::Pairs`](https://rdrr.io/pkg/S4Vectors/man/Pairs-class.html)
object with a copy of `x` in the `first` slot, and `x` with `step`
elements removed from its head and complemented with `NAs` in the
`second` slot.

## Details

This is a programming pattern that I use a lot in this package, to
prepare vectors fro parallel processing. I typically use it on metadata
columns of `GBreaks` objects. I created this function because I always
forget if I should use `head` or `tail`, which causes bugs or delays in
development.

## Examples

``` r
zipWithNext(LETTERS, step = 2)
#> Pairs object with 26 pairs and 0 metadata columns:
#>              first      second
#>        <character> <character>
#>    [1]           A           C
#>    [2]           B           D
#>    [3]           C           E
#>    [4]           D           F
#>    [5]           E           G
#>    ...         ...         ...
#>   [22]           V           X
#>   [23]           W           Y
#>   [24]           X           Z
#>   [25]           Y        <NA>
#>   [26]           Z        <NA>
```
