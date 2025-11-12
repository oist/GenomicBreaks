# Breakpoint proximity to Tandem Repeats

This function classifies breakpoints in terms of its proximity to tandem
repeats (but can be generalised to any property whose coverage is
binary).

## Usage

``` r
tan_bp(gr_ob, tan, tol, query_tf = FALSE)
```

## Arguments

- gr_ob:

  GRanges object containing pairwise alignment

- tan:

  GRanges object containing tandem repeat (or characteristic of
  interest) coverage

- tol:

  range in which breakpoints will be classified as "near" to a tandem
  repeat

- query_tf:

  logical value that, if set to TRUE, will classify bps in query
  coordinates as opposed to reference (default)

## Value

A list containing 3 GRanges objects of breakpoints far, near and within
tandem repeats.
