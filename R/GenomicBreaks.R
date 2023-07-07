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
#' @family Bioconductor API functions
#'
#' @importFrom methods callNextMethod new setClass
#'
#' @export GBreaks
#' @exportClass GBreaks

GBreaks <- setClass("GBreaks", contains = "GRanges")

setMethod("initialize", "GBreaks", function(.Object, ..., target = NULL, query = NULL, strand = NULL) {
  if (! (is.null(target) | is.null(query))) {
    gb <- target
    gb$query <- query
    if (! is.null(strand)) strand(gb) <- strand
    gb <- GBreaks(gb)
    return(gb)
  }
  gb <- callNextMethod(.Object, ...)
  if (length(gb) == 0) gb$query <- GRanges()
  gb
})

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
#'
#' @family Bioconductor API functions
#'
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

#' Test if a `GBreaks` object is sorted
#'
#' The only proper way to sort a [`GBreaks`] object is by ignoring strand, so
#' that inversions and deletions are easily detected and visualised.
#'
#' @param x A `GBreaks` object.
#'
#' @returns Returns `TRUE` or `FALSE`.
#'
#' @examples
#' isSorted(exampleInversion)
#'
#' @importFrom S4Vectors isSorted
#' @importMethodsFrom S4Vectors isSorted

setGeneric("isSorted", function(x) standardGeneric("isSorted"))

isSorted_GBreaks <- function(x)
  identical(x, sort(x, ignore.strand = TRUE))

#' @rdname isSorted
#' @export

setMethod("isSorted", "GBreaks", isSorted_GBreaks)

#' Range method for `GBreaks` objects
#'
#' This is a `range` method for [`GBreaks`] objects, that will run
#' [`GenomicRanges::range`] on its _target_ and _query_ ranges and will return
#' a new `GBreaks` object.
#'
#' @note `range` unconditionally ignores strand in `GBreaks` objects.
#'
#' @param x A `GBreaks` object.
#' @param with.revmap FALSE
#' @param ignore.strand FALSE
#' @param na.rm FALSE)
#' @param ... etc
#'
#' @returns tbd
#'
#' @family Bioconductor API functions
#'
#' @rdname range
#' @examples
#' range(exampleColinear3)
#' range(exampleDeletion)

range_GBreaks <- function(x, ..., with.revmap=FALSE, ignore.strand=FALSE, na.rm=FALSE) {
  gbl <- split(x, paste(seqnames(x), seqnames(x$query)), drop = TRUE)
  gb <- endoapply(gbl, \(gb) {
    rangeT <- range(granges(gb), ignore.strand = TRUE)
    rangeQ <- range(gb$query, ignore.strand = TRUE)
    GBreaks(target = rangeT, query = rangeQ)
  }) |> unlist()
  names(gb) <- NULL
  seqinfo(gb) <- seqinfo(x)
  seqinfo(gb$query) <- seqinfo(x$query)
  gb
}

#' @rdname range
#' @export

setMethod("range", "GBreaks", range_GBreaks)

#' subsetByOverlaps method for `GBreaks` objects
#'
#' This is a `subsetByOverlaps` method for [`GBreaks`] objects, that will run
#' [`GenomicRanges::subsetByOverlaps`] on its _target_ and _query_ ranges and
#' will return a new `GBreaks` object.
#'
#' @note `range` unconditionally ignores strand in `GBreaks` objects.
#'
#' @param x A `GBreaks` object.
#' @param ranges Another `GBreaks` object.
#' @param maxgap tbd
#' @param minoverlap tbd
#' @param type tbd
#' @param invert tbd
#' @param ... etc
#'
#' @returns tbd
#'
#' @family Bioconductor API functions
#'
#' @rdname subsetByOverlaps
#' @examples
#' subsetByOverlaps(exampleColinear3, exampleColinear3)

subsetByOverlaps_GBreaks <- function(x, ranges, maxgap=-1L, minoverlap=0L,
  type=c("any", "start", "end", "within", "equal"), invert=FALSE,  ...) {
  gb1 <- subsetByOverlaps(x, granges(ranges), maxgap = maxgap, minoverlap = minoverlap, type = type, invert = invert, ...)
  gb2 <- x[overlapsAny(x$query, ranges$query, minoverlap = minoverlap, type = type, ...)]
  sort(unique(c(gb1, gb2)))
}

#' @rdname subsetByOverlaps
#' @export
#'
setMethod("subsetByOverlaps", c("GBreaks", "GBreaks"), subsetByOverlaps_GBreaks)

