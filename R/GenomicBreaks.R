#' The GenomicBreaks class
#'
#' The `GBreaks` class is a simple wrapper to the [`GRanges`] class.
#'
#' @details Aligned sequences of the _target_ genome are represented as the main
#' ranges of the `GRanges` object, and their counterparts in the _query_ genome
#' are represented as a `GRanges` object sored in the the metadata column
#' `query`.
#'
#' @examples
#' gb <- load_genomic_breaks(system.file("extdata/SacCer3__SacPar.gff3.gz", package = "GenomicBreaks"))
#' gb
#' gb$query
#'
#' @importFrom methods new setClass
#'
#' @export

GBreaks <- setClass("GBreaks", contains = "GRanges")


#' Conversion from [`CNEr::Axt`] objects
#'
#' @note By default, in the _Axt_ objects produced by the [`CNEr::readAxt`]
#' function, the coordinates of the query genome are represented with the same
#' numerical value as in the original file.  In _Axt_ files, when the alignment
#' is on the minus strand, the query genome coordinates have their origin on
#' the reverse-complement strand.  Therefore, although in the _CNEr_ object
#' they are stored in a _GRanges_ object, they do not represent genomic ranges.
#' In order to do so, there is a [`CNEr::fixCoordinates`] function.
#'
#' @importFrom CNEr first second
#' @importFrom methods as setAs
#' @name as
#' @export

setAs("Axt", "GBreaks", function(from) {

  # First genome of axt object is target genome, second is query.
  gb       <- granges(first(from))
  gb$query <- granges(second(from))

  # In GBreaks object, strand information is carried by the target genome ranges
  strand(gb) <- strand(gb$query)
  strand(gb$query) <- "+"

  GBreaks(gb)
})
