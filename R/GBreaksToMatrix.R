#' Convert to binned matrix of hits
#'
#' @param gb A `GenomicBreaks` object.
#' @param ntile The number of bins in the matrix.
#'
#' @return A [`matrix`] object.
#'
#' @family Converter functions
#'
#' @examples
#' m <- GBreaksToMatrix(exampleColinear5, ntile = 10)
#' m
#' image(m)
#' makeOxfordPlots(exampleColinear5)
#'
#' @export

GBreaksToMatrix <- function(gb, ntile = 1000) {
  .guessSeqLengths <- function(gb) {
    zz <- guessSeqLengths(gb)
    if(is.vector(zz)) return(zz)
    if(!is.array(zz)) stop("unsupported class")
    zzz <- as.vector(zz)
    names(zzz) <- unlist(dimnames(zz))
    zzz
  }
  .tileGenome <- function(gb, ntile) {
    gb_seqLengths <- .guessSeqLengths(gb)
    gr <- tileGenome(gb_seqLengths, ntile = ntile - length(gb_seqLengths) + 1) |> unlist()
    names(gr) <- seq_along(gr)
    if (length(gr) != ntile) stop ("Could not generate the required number of tiles")
    gr
  }
  x <- .tileGenome(gb, ntile=ntile)
  #x <- c(x |> plyranges::mutate(strand = '+'), x |> plyranges::mutate(strand = '-')) |> sort (ignore.strand=T)
  y <- .tileGenome(gb$query, ntile=ntile)

  m <- matrix(0, nrow = ntile, ncol = ntile)

  for (i in 1:(ntile)) {
    bin <- x[i]
    o <- findOverlaps(gb, bin)
    oQuery <- findOverlaps(gb$query[queryHits(o)], y)
    m[i, subjectHits(oQuery)] <- 1
  }
  m
}

