#' Pairwise alignment of genomic ranges
#'
#' Retrieves DNA sequence from [`GBreaks`] or [`GRanges`] objects that are
#' properly related to a `BSgenome` package, and aligns them with the
#' [`pwalign::pairwiseAlignment()`] function.
#'
#' @returns Returns a [`pwalign::PairwiseAlignments`] object.
#'
#' @param pattern A `GBreaks` or a `GRanges` object
#' @param subject A `GBreaks` or a `GRanges` object
#' @param ... Additional arguments passed to `pairwiseAlignment`.
#'
#' @author Charles Plessy
#'
#' @examples
#' Scerevisiae <- BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae
#' # Very arbitrary example
#' gb <- GBreaks( target = GRanges("chrI: 1-20", seqinfo = seqinfo(Scerevisiae))
#'              , query  = GRanges("chrI:21-40", seqinfo = seqinfo(Scerevisiae)))
#' pairwiseAlignment(gb)
#'
#' @family Bioconductor API functions
#'
#' @importFrom pwalign pairwiseAlignment
#' @importMethodsFrom pwalign pairwiseAlignment

setGeneric("pairwiseAlignment", function(pattern, subject, ...) standardGeneric("pairwiseAlignment"))

pairwiseAlignment_GRanges_GRanges <- function(pattern, subject, ...) {
  seqP <- getSeq(pattern)
  seqS <- getSeq(subject)
  pairwiseAlignment(seqP, seqS, ...)
}

pairwiseAlignment_GBreaks <- function(pattern, subject, ...) {
  seqP <- getSeq(pattern)
  seqS <- getSeq(pattern$query)
  pairwiseAlignment(seqP, seqS, ...)
}

#' @rdname pairwiseAlignment
#' @export

setMethod("pairwiseAlignment", c("GRanges", "GRanges"), pairwiseAlignment_GRanges_GRanges)

#' @rdname pairwiseAlignment
#' @export

setMethod("pairwiseAlignment", c("GBreaks", NULL), pairwiseAlignment_GBreaks)
