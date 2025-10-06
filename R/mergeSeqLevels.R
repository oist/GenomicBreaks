#' Merge seqlevels in a larger one
#'
#' For scaffolding or plotting purposes, it may be useful to merge some
#' sequences into larger ones.
#'
#' @note Be careful that in some cases it is needed to "flip" the sequence
#' feature with [`reverse`] before merging, for instance when colinearity is
#' with its reverse strand.
#'
#' @param gr A [`GenomicRanges::GRanges`] object.
#' @param seqs A character vector of [`Seqinfo::seqlevels`] from `gr`
#' @param name The name of the new sequence level to be added
#'
#' @returns Returns a modified `GRanges` object in which the sequences have been
#' merged.  Its [`Seqinfo::seqinfo`] has a new entry for the new level, and the old
#' levels are not removed.  If no `seqlengths` were present in the original
#' object, they are arbitrarily set as the maximal end value for each `seqlevel`.
#'
#' @returns The `mergeSeqLevels_to_DF` function returns a `DataFrame` in which
#' the `start` and `end` columns are in `numeric` mode.  This is to cirvumvent
#' the fact that `GenomicRanges` object hardcode the mode of _start_ and _end_
#' positions to `integer`, which does not allow values larger than
#' 2,147,483,647, which does not allow to merge sequence levels of mammalian
#' or larger-scale genomes.
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
#' @family scaffolding functions
#'
#' @author Charles Plessy
#'
#' @include guessSeqLenghts.R
#'
#' @export

mergeSeqLevels <- function(gr, seqs, name) {
  if (identical(seqs, name)) return(gr)

  # Calculate how much to add to the coordinates of each seqlevel before merging
  if (all(is.na(seqlengths(gr))))
    gr <- forceSeqLengths(gr)
  lengths <- seqlengths(gr)[seqs]
  mode(lengths) <- "numeric"                      # 32-bit integers are too small for
  addlengths <- c(cumsum(head(lengths, -1)))      # the cumulative sum on this line
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

#' @rdname mergeSeqLevels
#' @export

mergeSeqLevels_to_DF <- function(gr, seqs, name) {
  # Calculate how much to add to the coordinates of each seqlevel before merging
  if (all(is.na(seqlengths(gr))))
    gr <- forceSeqLengths(gr)
  lengths <- seqlengths(gr)[seqs]
  mode(lengths) <- "numeric"                      # 32-bit integers are too small for
  addlengths <- c(cumsum(head(lengths, -1)))      # the cumulative sum on this line
  names(addlengths) <- tail(names(lengths), -1)

  # A DataFrame representing the GRanges, but using numeric mode instead of integer
  DF <- DataFrame(
    seqnames = seqnames(gr),
    start = start(gr) |> as.numeric(),
    end = end(gr) |> as.numeric(),
    strand = strand(gr))

  levels(DF$seqnames) <- c(levels(DF$seqnames), name)

  # Increment coordinates of each contigs and rename them
  for(contig in names(addlengths)) {
    DF[DF$seqnames == contig, "end"]      <- DF[DF$seqnames == contig, "end"]   + addlengths[contig]
    DF[DF$seqnames == contig, "start"]    <- DF[DF$seqnames == contig, "start"] + addlengths[contig]
    DF[DF$seqnames == contig, "seqnames"] <- name
  }
  DF[DF$seqnames == seqs[1], "seqnames"] <- name

  # Return the DataFrame
  DF
}
