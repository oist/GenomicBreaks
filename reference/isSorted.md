# Test if a `GBreaks` object is sorted

The only proper way to sort a
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
object is by ignoring strand, so that inversions and deletions are
easily detected and visualised.

## Usage

``` r
isSorted(x)

# S4 method for class 'GBreaks'
isSorted(x)
```

## Arguments

- x:

  A `GBreaks` object.

## Value

Returns `TRUE` or `FALSE`.

## Examples

``` r
isSorted(exampleInversion)
#> [1] TRUE
```
