# Load pairwise genome alignments

Loads alignments of a *query* genome to a *target* genome from a text
file in *General Feature Format 3* (GFF3) or *Multiple Alignemnt Format*
(MAF). By convention, the *target* genome is the one that was indexed by
the aligner.

## Usage

``` r
load_genomic_breaks(
  file,
  target_bsgenome = NULL,
  query_bsgenome = NULL,
  sort = TRUE,
  type = "match_part"
)
```

## Arguments

- file:

  Path to a file in GFF3 or MAF format. The file can be compressed with
  *gzip*.

- target_bsgenome:

  A `BSgenome` object or path to a FASTA file representing the *target*
  genome.

- query_bsgenome:

  A `BSgenome` object or path to a FASTA file representing the *query*
  genome.

- sort:

  Returns the object sorted, ignoring strand information.

- type:

  In GFF3 files, *Sequence Ontology* term representing an alignment
  block (default: `match_part`).

## Value

Returns a
[`GBreaks`](https://oist.github.io/GenomicBreaks/reference/GBreaks-class.md)
object where each element represents a pairwise alignment block. The
`granges` part of the object contains the coordinates on the *target*
genome, and the `query` metadata column contains the *query* coordinates
in `GRanges` format. The `seqinfo` of each `BSgenome` object passed as
parameters are copied to the `GRanges` objects accordingly.

## Details

When the input is in GFF3 files, this function expects the pairwise
alignment to be represented in in the following way:

- Alignments blocks are represented by entries in specific sequence
  ontology term in the *type* column. Other entries will be discarded.
  The default type is `match_part`.

- The coordinate system of the file is the one of the *target* genome.

- The `Target` tag in the *attribute* column contains the coordinates of
  the match in the *query* genome. (Sorry that it is confusing…)

- Stand information is set so that *query* genome coordinates are always
  on the *plus* strand.

## See also

The [MAF format](http://www.genome.ucsc.edu/FAQ/FAQformat.html#format5)
documentation on the UCSC genome browser website, and the [GFF3
specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
from the *Sequence Ontology* group.

Other Data loading functions:
[`readMAF()`](https://oist.github.io/GenomicBreaks/reference/readMAF.md),
[`readTrainFile()`](https://oist.github.io/GenomicBreaks/reference/readTrainFile.md)

## Examples

``` r
load_genomic_breaks(system.file("extdata/contigs.genome.maf.gz", package = "GenomicBreaks"))
#> GBreaks object with 2 ranges and 4 metadata columns:
#>         seqnames     ranges strand |     score
#>            <Rle>  <IRanges>  <Rle> | <integer>
#>   [1] MT192765.1    25-8666      + |     52990
#>   [2] MT192765.1 8882-29829      - |    128566
#>                                             query   aLength   matches
#>                                         <GRanges> <integer> <integer>
#>   [1]    NODE_2_length_8774_cov_178.827802:6-8647      8642      8642
#>   [2] NODE_1_length_20973_cov_191.628754:26-20973     20948     20945
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome

if (FALSE) { # \dontrun{
library("BSgenome.Scerevisiae.UCSC.sacCer3")
load_genomic_breaks(
  system.file("extdata/SacCer3__SacPar.gff3.gz", package = "GenomicBreaks"),
  target = Scerevisiae,
  query = NULL)
} # }
```
