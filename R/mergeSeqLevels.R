#' Merge seqlevels in a larger one
#'
#' For scaffolding or plotting purposes, it may be useful to merge some
#' sequences into larger ones.
#'
#' @note Be careful that in some cases it is needed to "flip" the sequence
#' feature before merging, if it maps colinearly with its reverse strand.
#'
#' @param gr A `GRanges` object.
#' @param seqs A character vector of [`seqlevels`] from `gr`
#' @param name The name of the new sequence level to be added
#'
#' @returns Returns a modified `GRanges` object in which the sequences have been
#' merged.  Its [`seqinfo`] has a new entry for the new level, and the old
#' levels are not removed.
#'
#' @examples
#' gb       <- GRanges(c("XSR:101-180:+", "XSR:201-300:+",  "XSR:320-400:+"))
#' gb$query <- GRanges(c( "S1:101-200",      "S2:1-100",    "S3:1-100"))
#' seqlengths(gb$query) <- c(200, 100, 100)
#' genome(gb$query) <- "GenomeX"
#' isCircular(gb$query) <- rep(FALSE, 3)
#' seqinfo(gb$query)
#' gb <- GBreaks(gb)
#' gb$query <- mergeSeqLevels(gb$query, c("S2", "S3"), "Scaf1")
#' gb
#' seqinfo(gb$query)
#'
#' mergeSeqLevels(gb, seqlevelsInUse(gb), "AllMerged")
#'
#' @family modifier functions
#'
#' @author Charles Plessy
#'
#' @export

mergeSeqLevels <- function(gr, seqs, name) {

  # Calculate how much to add to the coordinates of each seqlevel before merging
  lengths <- seqlengths(seqinfo(gr))[seqs]
  addlengths <- c(cumsum(head(lengths, -1)))
  names(addlengths) <- tail(names(lengths), -1)

  # Add a new seqlevel
  genomeName <- unique(genome(gr))
  seqlevels(gr) <- c(name, seqlevels(gr))
  seqlengths(seqinfo(gr))[name] <- sum(lengths)
  isCircular(seqinfo(gr))[name] <- FALSE
  genome(seqinfo(gr))[name] <- genomeName

  # Increment coordinates of each contigs and rename them
  for(contig in names(addlengths)) {
    end     (gr[seqnames(gr) == contig]) <- end    (gr[seqnames(gr) == contig]) + addlengths[contig]
    start   (gr[seqnames(gr) == contig]) <- start  (gr[seqnames(gr) == contig]) + addlengths[contig]
    seqnames(gr[seqnames(gr) == contig]) <- factor(name, levels = seqlevels(gr))
  } |> suppressWarnings()
  seqnames(gr[seqnames(gr) == seqs[1]]) <- factor(name, levels = seqlevels(gr))

  # REturn the modified object
  gr
}
