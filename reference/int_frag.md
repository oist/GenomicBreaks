# Integrated intersection of both query and reference subsets on a continuous alignment

Integrated intersection of both query and reference subsets on a
continuous alignment

## Usage

``` r
int_frag(og_r, og_q, sub_r, sub_q)
```

## Arguments

- og_r:

  GRanges object containing the refernce genome information of the
  original alignment

- og_q:

  GRanges object containing the query genome information of the original
  alignment

- sub_r:

  GRanges object containing the reference subset

- sub_q:

  GRanges object containing the query subset

## Value

GRanges object of alignment intersected on both reference and query
side. Query data held in metadata column
