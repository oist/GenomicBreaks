# Read a MAF file

Reads a pairwise genome alignment in MAF format. The file can be plain
text or compressed with `gzip`.

## Usage

``` r
readMAF(inputFileName)
```

## Arguments

- inputFileName:

  The name of the file to read

## Value

a `list` object with containg coordinates of the alignments in both
genomes and other information such as alignment width and number of
matches.

## Details

Known limitations: Does not expand shell metacharacters. Trusts blindly
file extension to determine compression. Does not perform any validation
on the file format. Assumes that the score comes first in the 'a' lines.

## See also

Other Data loading functions:
[`load_genomic_breaks()`](https://oist.github.io/GenomicBreaks/reference/load_genomic_breaks.md),
[`readTrainFile()`](https://oist.github.io/GenomicBreaks/reference/readTrainFile.md)
