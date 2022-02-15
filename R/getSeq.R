#' Get genomic sequences
#'
#' Extract the sequences of the genomic ranges from a [`GenomicRanges::GRanges`]
#' or a [`GBreaks`] object.
#'
#' @returns A [`Biostrings::DNAStringSet`] object containing the extracted
#' sequence(s).
#'
#' @param x A `GBreaks` or a `GRanges` object.
#' @param ... Extra arguments (not used, but present for compatibility).
#'
#' @author Charles Plessy
#'
#' @examples
#' Scerevisiae <- BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae
#' getSeq(GRanges("chrI:1-20", seqinfo = seqinfo(Scerevisiae)))
#'
#' @family Bioconductor API functions
#' @seealso See also the [`Biostrings::getSeq`] function.
#'
#' @importFrom Biostrings getSeq
#' @importFrom BSgenome getBSgenome
#' @importMethodsFrom Biostrings getSeq

setGeneric("getSeq", function(x, ...) standardGeneric("getSeq"))

getSeq_GRanges <- function(x, ...) {
  genome <- unique(genome(x))
  if (any(is.na(genome)))
    stop("Can not getSeq objets where `genome` is not defined.")
  getSeq(getBSgenome(genome), x)
}

getSeq_GBreaks <- function(x, ...) {
  getSeq(granges(x))
}

#' @rdname getSeq
#' @export

setMethod("getSeq", "GRanges", getSeq_GRanges)

#' @rdname getSeq
#' @export

setMethod("getSeq", "GBreaks", getSeq_GBreaks)
