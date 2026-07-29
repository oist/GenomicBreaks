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
#' @param target_bsgenome A `BSgenome` object or path to a FASTA file representing the _target_ genome.
#' @param query_bsgenome A `BSgenome` object or path to a FASTA file representing the _query_ genome.
#' @param sort Returns the object sorted, ignoring strand information.
#' @param type In GFF3 files, _Sequence Ontology_ term representing an alignment
#'        block (default: `match_part`).
#'
#' @family Data loading functions
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
#' load_genomic_breaks(system.file("extdata/contigs.genome.maf.gz", package = "GenomicBreaks"))
#' load_genomic_breaks(system.file("extdata/MT192765___genome-OY074094.blasttabplus.gz",
#'   package = "GenomicBreaks"))
#'
#' \dontrun{
#' library("BSgenome.Scerevisiae.UCSC.sacCer3")
#' load_genomic_breaks(
#'   system.file("extdata/SacCer3__SacPar.gff3.gz", package = "GenomicBreaks"),
#'   target = Scerevisiae,
#'   query = NULL)
#' }
#'
#' @importFrom utils read.delim
#' @importFrom rtracklayer import.gff3
#' @importFrom Biostrings readDNAStringSet
#' @export

load_genomic_breaks <- function (
    file,
    target_bsgenome = NULL,
    query_bsgenome = NULL,
    sort = TRUE,
    type = "match_part")
{
  file <- normalizePath(file, mustWork = TRUE)
  if (grepl(".gff3$|gff3.gz$|.gff$|gff.gz$", file)) {
    load_genomic_breaks_function <- load_genomic_breaks_GFF
  } else if (grepl(".maf$|.maf.gz$", file)) {
    load_genomic_breaks_function <- load_genomic_breaks_MAF
  } else if (grepl("blasttabplus$|blasttabplus.gz$", file)) {
    load_genomic_breaks_function <- load_genomic_breaks_blasttabplus
  } else {

    stop("Unsupported file type: extension should be gff, gff.gz, gff3, gff3.gz, maf, or maf.gz.")
  }
  gb <- load_genomic_breaks_function(
    file = file, target_bsgenome = target_bsgenome,
    query_bsgenome = query_bsgenome, sort = sort, type = type)
  gb
}

load_genomic_breaks_blasttabplus <- function (
    file,
    target_bsgenome = NULL,
    query_bsgenome = NULL,
    sort = TRUE,
    type = "match_part")
{
  blasttabplusformat <- rbind(
    c(name = "qname", type = "factor"),
    c(       "tname",        "factor"),
    c(       "P",            "numeric"),
    c(       "alength",      "integer"),
    c(       "mismatches",   "integer"),
    c(       "gaps",         "integer"),
    c(       "qstart",       "integer"),
    c(       "qend",         "integer"),
    c(       "tstart",       "integer"),
    c(       "tend",         "integer"),
    c(       "E",            "numeric"),
    c(       "bit",          "numeric"),
    c(       "qlen",         "integer"),
    c(       "tlen",         "integer"),
    c(       "score",        "integer")
  ) |> as.data.frame()
  df <- read.delim( file
                  , header = FALSE
                  , colClasses = blasttabplusformat$type
                  , col.names = blasttabplusformat$name)

  toSeqinfo <- function(si)
    DataFrame( seqlengths = si[,2]
               , row.names  = si[,1]
               , isCircular = FALSE
               , genome = NA_character_) |> as("Seqinfo")

  mkSeqinfo <- function(genome, df)
    switch( class(genome)
          , "NULL" = toSeqinfo(unique(df))
          , "character" = seqinfo(readDNAStringSet(genome))
          , seqinfo(genome)
    )

  mkGranges <- function(seqnames, start, end, seqinfo = NULL) {
    flip <- end < start
    GRanges( seqnames = seqnames
           , ranges   = IRanges( start = ifelse(flip, end, start)
                               , end   = ifelse(flip, start, end))
           , strand   = ifelse(flip, "-", "+")
           , seqinfo  = seqinfo
    )
  }

  grt <- mkGranges(df$tname, df$tstart, df$tend, mkSeqinfo(target_bsgenome, df[,c("tname", "tlen")]))
  grq <- mkGranges(df$qname, df$qstart, df$qend, mkSeqinfo(query_bsgenome,  df[,c("qname", "qlen")]))

  stopifnot(all(strand(grq)=="+"))
  strand(grq) <- "*"

  # Let's not report P as it is lossy and can not be used to recover number of mismatches.
  #
  # Quoting from the maf-convert source code:
  #
  # alignmentColumns = list(zip(rowA, rowB))
  # alnSize = len(alignmentColumns)
  # matches = sum(x == y for x, y in alignmentColumns)
  # mismatches = alnSize - matches - rowA.count("-") - rowB.count("-")
  # matchPercent = "%.2f" % (100.0 * matches / alnSize)

  gb <- GBreaks (target = grt, query = grq)
  score(gb)     <- df$score
  gb$alength    <- df$alength
  gb$mismatches <- df$mismatches
  gb$gaps       <- df$gaps

  if (sort) gb <- sort(gb, ignore.strand = TRUE)
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
  if (!is.null(target_bsgenome)) {
    if (is.character(target_bsgenome)) {
      target_seqinfo <- seqinfo(readDNAStringSet(target_bsgenome))
    } else {
      target_seqinfo <- seqinfo(target_bsgenome)
    }
    gb <- GRanges(gb, seqinfo = target_seqinfo)
  }
  # Discard cross_genome_match parent (used for block display in Zenbu)
  gb <- gb[gb$type == type]
  # Discard unused information
  gb$phase <- gb$Parent <- gb$Target <- gb$ID <- gb$source <- gb$type <- NULL
  # Convert query coordinates to GRanges
  if (!is.null(query_bsgenome)) {
    if (is.character(query_bsgenome)) {
      query_seqinfo <- seqinfo(readDNAStringSet(query_bsgenome))
    } else {
      query_seqinfo <- seqinfo(query_bsgenome)
    }
    gb$query <- GRanges(gb$Name, seqinfo = query_seqinfo)
  } else {
    gb$query <- GRanges(gb$Name)
  }
  gb$Name <- NULL
  if (sort) gb <- sort(gb, ignore.strand = TRUE)
  as(gb, "GBreaks")
}

load_genomic_breaks_MAF <- function (
    file,
    target_bsgenome = NULL,
    query_bsgenome = NULL,
    sort = TRUE,
    type = NULL)
{
  l <- readMAF(file)
  # Convert coordinates to GenomicRanges system
  # http://www.genome.ucsc.edu/FAQ/FAQformat.html#format5
  # Make minus-strand coordinates relative to the start of the sequence
  l$start1 <- ifelse(l$strand1 == "+",
                     l$start1,
                     l$seqlengths1 - l$start1 - l$length1)
  l$start2 <- ifelse(l$strand2 == "+",
                     l$start2,
                     l$seqlengths2 - l$start2 - l$length2)
  # Add 1 to the starts
  # http://genomewiki.ucsc.edu/index.php/Coordinate_Transforms
  l$start1 = l$start1 + 1;
  l$start2 = l$start2 + 1;
  # Build GBreaks object
  # Note that GRanges(...,strand=FALSE) gives +, and strand=TRUE gives -.
  gb <- GRanges(l$seqnames1, IRanges(l$start1, width = l$length1), strand = l$strand1 != l$strand2)
  score(gb) <- l$scores
  seqlengths1        <- l $ seqlengths1 [!duplicated(l$seqnames1)]
  names(seqlengths1) <- l $ seqnames1   [!duplicated(l$seqnames1)]
  seqlengths(gb) <- seqlengths1
  gb$query <- GRanges(l$seqnames2, IRanges(l$start2, width = l$length2))
  seqlengths2        <- l $ seqlengths2 [!duplicated(l$seqnames2)]
  names(seqlengths2) <- l $ seqnames2   [!duplicated(l$seqnames2)]
  seqlengths(gb$query) <- seqlengths2
  gb$aLength <- l$aLength
  gb$matches <- l$matches
  if (sort) gb <- sort(gb, ignore.strand = TRUE)
  as(gb, "GBreaks")
}
