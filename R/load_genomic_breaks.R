#' Load pairwise genome alignments
#'
#' Loads alignments of a _query_ genome to a _target_ genome from a text file in
#' GFF3 format.  By convention, the _target_ genome is the one that was indexed
#' by the aligner.
#'
#' This function expects the pairwise alignment to be represented in GFF3 format
#' in the following way:
#'
#' * Alignments blocks are represented by entries in specific sequence ontology
#'   term in the _type_ column.  Other entries will be discarded.  The default
#'   type is `match_part`.
#' * The coordinate system of the file is the one of the _target_ genome.
#' * The `Target` tag in the _attribute_ column contains the coordinates of
#'   the match in the _query_ genome.  (Sorry that it is confusing…)
#' * Stand information is set so that _query_ genome coordinates are always on
#'   the _plus_ strand.
#'
#' @param file Path to a file in GFF3 format.
#' @param target_bsgenome A `BSgenome` object representing the target genome.
#' @param query_bsgenome A `BSgenome` object representing the query genome.
#' @param sort Returns the object sorted, ignoring strand information.
#' @param type Sequence ontology term representing an alignment block
#'        (default: `match_part`).
#'
#' @return Returns a [`GRanges`] object where each element represents a pairwise
#' alignment block.  The `granges` part of the object contains the coordinates
#' on the _target_ genome, and the `query` metadata column contains the _query_
#' coordinates in `GRanges` format.  The `seqinfo` of each `BSgenome` object
#' passed as parameters are copied to the `GRanges` objects accordingly.
#'
#' @examples
#' \dontrun{
#' library("BSgenome.Odioica.local.OSKA2016")
#' library("BSgenome.Odioica.local.OKI2018.I69")
#' load_genomic_breaks(
#'   system.file("extdata/OSKA2016__I69-5.gff3.gz", package="GenomicBreaks"),
#'   target = OSKA2016,
#'   query = OKI2018_I69)
#' }
#'
#' @importFrom rtracklayer import.gff3
#' @export

load_genomic_breaks <- function (
  file,
  target_bsgenome = NULL,
  query_bsgenome = NULL,
  sort = TRUE,
  type = "match_part")
{
  gb <- import.gff3(file)
  if (! is.null(target_bsgenome))
    gb <- GRanges(gb, seqinfo = seqinfo(target_bsgenome))
  # Discard cross_genome_match parent (used for block display in Zenbu)
  gb <- gb[gb$type == type]
  # Discard unused information
  gb$phase <- gb$Parent <- gb$Target <- gb$ID <- gb$source <- gb$type <- NULL
  # Convert query coordinates to GRanges
  if (! is.null(query_bsgenome))
    gb$query <- GRanges(gb$Name, seqinfo = seqinfo(query_bsgenome))
  else
    gb$query <- GRanges(gb$Name)
  gb$Name <- NULL
  if (sort) gb <- sort(gb, ignore.strand = TRUE)
  gb
}
