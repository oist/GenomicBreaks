#' Load pairwise genome alignments
#'
#' Loads alignments of a _query_ genome to a _target_ genome from a text file in
#' _General Feature Format 3_ (GFF3) or _Multiple Alignemnt Format_ (MAF).  By
#' convention, the _target_ genome is the one that was indexed by the aligner.
#'
#' When the input is in GFF3 files, this function expects the pairwise alignment
#' to be represented in in the following way:
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
#' @param file Path to a file in GFF3 or MAF format.  The file can be compressed
#'        with _gzip_.
#' @param target_bsgenome A `BSgenome` object representing the _target_ genome.
#' @param query_bsgenome A `BSgenome` object representing the _query_ genome.
#' @param sort Returns the object sorted, ignoring strand information.
#' @param type In GFF3 files, _Sequence Ontology_ term representing an alignment
#'        block (default: `match_part`).
#'
#' @seealso The [MAF format](http://www.genome.ucsc.edu/FAQ/FAQformat.html#format5)
#' documentation on the UCSC genome browser website, and the
#' [GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
#' from the _Sequence Ontology_ group.
#'
#' @return Returns a [`GBreaks`] object where each element represents a pairwise
#' alignment block.  The `granges` part of the object contains the coordinates
#' on the _target_ genome, and the `query` metadata column contains the _query_
#' coordinates in `GRanges` format.  The `seqinfo` of each `BSgenome` object
#' passed as parameters are copied to the `GRanges` objects accordingly.
#'
#' @examples
#' \dontrun{
#' library("BSgenome.Scerevisiae.UCSC.sacCer3")
#' load_genomic_breaks(
#'   system.file("extdata/SacCer3__SacPar.gff3.gz", package = "GenomicBreaks"),
#'   target = Scerevisiae,
#'   query = NULL)
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
  file <- normalizePath(file, mustWork = TRUE)
  if (grepl(".gff3$|gff3.gz$", file)) {
    load_genomic_breaks_function <- load_genomic_breaks_GFF
  } else if (grepl(".maf$|.maf.gz$", file)) {
    load_genomic_breaks_function <- load_genomic_breaks_MAF
  }
  gb <- load_genomic_breaks_function(
    file = file, target_bsgenome = target_bsgenome,
    query_bsgenome = query_bsgenome, sort = sort, type = type)
  gb
}

load_genomic_breaks_GFF <- function (
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
  as(gb, "GBreaks")
}

load_genomic_breaks_MAF <- function (
    file,
    target_genome = NULL,
    query_genome = NULL,
    sort = TRUE,
    type = NULL)
{
  l <- readMAF(file)
  # Convert coordinates to GenomicRanges system
  # http://www.genome.ucsc.edu/FAQ/FAQformat.html#format5
  # Make minus-strand coordinates relative to the start of the sequence
  l$start2 <- ifelse(l$strand == "+",
                     l$start2,
                     l$seqlengths2 - l$start2 - l$length2 + 1)
  # Add 1 to the starts
  # http://genomewiki.ucsc.edu/index.php/Coordinate_Transforms
  l$start1 = l$start1 + 1;
  l$start2 = l$start2 + 1;
  # Build GBreaks object
  gb <- GRanges(l$seqnames1, IRanges(l$start1, width = l$length1), strand = l$strand)
  score(gb) <- l$scores
  seqlengths1 <- unique(l$seqlengths1)
  names(seqlengths1) <- unique(l$seqnames1)
  seqlengths(gb) <- seqlengths1
  gb$query <- GRanges(l$seqnames2, IRanges(l$start2, width = l$length2))
  seqlengths2 <- unique(l$seqlengths2)
  names(seqlengths2) <- unique(l$seqnames2)
  seqlengths(gb$query) <- seqlengths2
  gb$aLength <- l$aLength
  gb$matches <- l$matches
  if (sort) gb <- sort(gb, ignore.strand = TRUE)
  as(gb, "GBreaks")
}
