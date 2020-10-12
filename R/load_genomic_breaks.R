#' Load pairwise genome alignments
#'
#' Loads alignments of a query genome to a reference genome from a text file in
#' GFF3 format
#'
#' This function expects the pairwise alignment to be represented in GFF3 format
#' in the following way:
#'
#' * Alignments blocks are represented by `match_part` entries in the _type_
#'   column  Other entries will be discarded.
#' * The coordinate system of the file is the one of the _reference_ genome.
#' * The `Target` tag in the _attribute_ column contains the coordinates of
#'   the match in the _query_ genome.
#' * Stand information is set so that _query_ genome coordinates are always on
#'   the _plus_ strand.
#'
#' @param file Path to a file in GFF3 format
#' @param bsgenome_ref A `BSgenome` object representing the reference genome.
#' @param bsgenome_query A `BSgenome` object representing the query genome.
#'
#' @return Returns a `GRanges object` where each element represents a pairwise
#' alignment block.  The `granges` part of the object contains the coordinates
#' on the _reference_ genome, and the `query` metadata column contains the query
#' coordinates in `GRanges` format.  The `seqinfo` of each `BSgenome` object
#' passed as parameters are copied to the `GRanges` objects accordingly.
#'
#' @importFrom rtracklayer import.gff3
#' @export

load_genomic_breaks <- function(file, bsgenome_ref, bsgenome_query) {
  gb <- import.gff3(file, genome = deparse(substitute(bsgenome_ref)))
  # Discard cross_genome_match parent (used for block display in Zenbu)
  gb <- gb[gb$type == "match_part"]
  # Discard unused information
  gb$phase <- gb$Parent <- gb$Target <- gb$ID <- gb$source <- gb$type <- NULL
  # Convert query coordinates to GRanges
  gb$query <- GRanges(gb$Name, seqinfo = seqinfo(bsgenome_query))
  gb$Name <- NULL
  gb
}
