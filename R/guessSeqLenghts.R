#' Guesstimate seqlevel lengths
#'
#' When no [`GenomeInfoDb::seqlengths`] are available, one can resort to set
#' them as to the maximal end coordinate found in the object.
#'
#' @param gr A [`GenomicRanges::GRanges`] object
#'
#' @returns Returns sequence lengths that have been guessed as the maximal end
#' coordinates found in the `gr`, or the `gr`'s `seqlengths` if if they did
#' already exist.
#'
#' @author Charles Plessy
#'
#' @family modifier functions
#'
#' @examples
#' exampleTranslocation
#' seqlengths(exampleTranslocation$query)
#' guessSeqLengths(exampleTranslocation$query)
#' gb2 <- exampleTranslocation
#' seqlengths(gb2$query) <- NA
#' guessSeqLengths(gb2$query)
#'
#' @export

guessSeqLengths <- function(gr) {
  if (!any(is.na(seqlengths(gr)))) return (seqlengths(gr))
  tapply(end(gr), seqnames(gr), max)
}

#' Set sequence lengths if there were none
#'
#' Takes a [`GBreaks`] or a [`GenomicRanges::GRanges`] object and if `seqlengths` are not
#' available, sets them using the output of the [`guessSeqLengths`] function.
#'
#' @param x A `GBreaks` or a `GRanges` object.
#'
#' @returns [`forceSeqLengths`] returns the object in which the sequence lengths
#' have been set to the maximal end coordinates found in the object if if they
#' did not exist.  For `GBreaks` objects it handles both the _target_ and the
#' _query_ ranges at the same time.
#'
#' @family modifier functions
#'
#' @examples
#' # Prepare an example object with no seqlengths
#' gb <- exampleTranslocation
#' seqlengths(gb) <-  seqlengths(gb$query) <- NA
#' gb
#'
#' # Note that the new seqlengths returned by forceSeqLengths are shorter
#' # because we can not guess about length of the unaligned ends.
#' forceSeqLengths(gb)  |> seqlengths()
#' exampleTranslocation |> seqlengths()
#'
#' # forceSeqLengths can take whole GBreaks objects as input, or simple GRanges
#' forceSeqLengths(gb)$query    |> seqlengths()
#' forceSeqLengths(gb$query)    |> seqlengths()
#' forceSeqLengths(granges(gb)) |> seqlengths()
#'
#' @export

setGeneric("forceSeqLengths", function(x) standardGeneric("forceSeqLengths"))

forceSeqLengths_GRanges <- function(x) {
  if (!any(is.na(seqlengths(x)))) return (x)
  seqlengths(x) <- as.vector(guessSeqLengths(x))
  x
}

forceSeqLengths_GBreaks <- function(x) {
  if (any(is.na(seqlengths(x))))
    seqlengths(x) <- as.vector(guessSeqLengths(x))
  if (any(is.na(seqlengths(x$query))))
    seqlengths(x$query) <- as.vector(guessSeqLengths(x$query))
  x
}

#' @rdname forceSeqLengths

setMethod("forceSeqLengths", "GRanges", forceSeqLengths_GRanges)

#' @rdname  forceSeqLengths


setMethod("forceSeqLengths", "GBreaks", forceSeqLengths_GBreaks)
