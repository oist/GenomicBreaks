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
load_genomic_breaks(system.file("extdata/SacCer3__SacPar.blasttabplus.gz",
  package = "GenomicBreaks"))
#> GBreaks object with 677 ranges and 6 metadata columns:
#>         seqnames        ranges strand |                       query     score
#>            <Rle>     <IRanges>  <Rle> |                   <GRanges> <integer>
#>     [1]     chrI     2452-2714      - |   NC_047500.1:770507-770785       736
#>     [2]     chrI     6952-9750      + |      NC_047487.1:9556-12401      9367
#>     [3]     chrI    9779-10007      + |     NC_047487.1:12486-12713       216
#>     [4]     chrI   11455-11624      + |     NC_047487.1:13352-13522       549
#>     [5]     chrI   11625-12502      - |   NC_047487.1:202778-203667      2937
#>     ...      ...           ...    ... .                         ...       ...
#>   [673]   chrXVI 921147-928804      + |   NC_047502.1:881568-889175     23396
#>   [674]   chrXVI 928893-929583      + |   NC_047502.1:889416-890120       619
#>   [675]   chrXVI 929848-933046      + |   NC_047502.1:890289-893554      7363
#>   [676]   chrXVI 933405-936088      + |   NC_047502.1:893560-896238     10002
#>   [677]   chrXVI 936726-942083      + | NC_047493.1:1107170-1112502     25819
#>                 P   alength mismatches      gaps
#>         <numeric> <integer>  <integer> <integer>
#>     [1]     80.65       279         38         5
#>     [2]     84.45      2856        377        22
#>     [3]     68.26       230         70         3
#>     [4]     84.88       172         23         2
#>     [5]     79.57       920        116         3
#>     ...       ...       ...        ...       ...
#>   [673]     82.12      7769       1117        82
#>   [674]     66.94       738        164        22
#>   [675]     74.89      3345        615        39
#>   [676]     87.50      2704        293        19
#>   [677]     95.53      5370        191         6
#>   -------
#>   seqinfo: 17 sequences from an unspecified genome

if (FALSE) { # \dontrun{
library("BSgenome.Scerevisiae.UCSC.sacCer3")
load_genomic_breaks(
  system.file("extdata/SacCer3__SacPar.gff3.gz", package = "GenomicBreaks"),
  target = Scerevisiae,
  query = NULL)
} # }
```
