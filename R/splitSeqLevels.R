#' Split a seqlevel in two pieces
#'
#' For scaffolding or plotting purposes, it may be useful to split some
#' sequences into smaller ones.
#'
#' @note This function only splits at breakpoints.
#'
#' @param gb A `GBreaks` object.
#' @param seq A one seqlevel from `gr`
#' @param bp The breakpoint where to split.
#'
#' @returns Returns a modified `GRanges` object in which the sequence has been
#' split  Its [`GenomeInfoDb::seqinfo`] has a new entry for the new levels, and the old
#' level is not removed.  If no `seqlengths` were present in the original
#' object, they are arbitrarily set as the maximal end value for each `seqlevel`.
#'
#' @examples
#' splitSeqLevel(exampleInsertion, "chrA", 200)
#'
#' @family modifier functions
#' @family scaffolding functions
#'
#' @author Charles Plessy
#'
#' @include guessSeqLenghts.R
#'
#' @export

splitSeqLevel <- function(gb, seq, bp) {
  # Isolate the target ranges that are in the selected seqlevel.
  gbl.all <- split(gb, seqnames(gb))
  gb <- gbl.all[[seq]]

  # Prepare new levels and lengths
  seq_1 <- paste0(seq, "_1")
  seq_2 <- paste0(seq, "_2")
  seqlevels(gb) <- c(seqlevels(gb), seq_1, seq_2)
  seqlengths(gb)[seq_1] <- bp
  seq_length <- seqlengths(gb)[[seq]]
  if(is.na(seq_length)) seq_length <- guessSeqLengths(gb)[[seq]]
  seqlengths(gb)[seq_2] <- seq_length - bp
  genome(gb)[c(seq_1, seq_2)] <- unique(genome(gbl.all))

  # Split, change, merge
  gbl <- split(gb, ifelse(end(gb) <= bp, seq_1, seq_2))
  seqnames(gbl[[seq_1]]) <- factor(seq_1, levels = seqlevels(gb))
  gbl[[seq_2]] <- shift(gbl[[seq_2]], - bp)
  seqnames(gbl[[seq_2]]) <- factor(seq_2, levels = seqlevels(gb))
  gbl.all[[seq]] <- unlist(gbl)

  # Return new GBreaks object
  sort(ignore.strand = TRUE, unlist(gbl.all))
}
