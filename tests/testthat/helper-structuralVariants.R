# Move the first or the last range of a GBreaks object to a new target
# sequence, ordered before or after the existing ones so that the object stays
# sorted.  The query ranges are unchanged.
moveToOtherTargetSeq <- function(gb, i) {
  stopifnot(i %in% c(1, length(gb)))
  if (i == 1) {
    seqlevels(gb) <- c("chr0", seqlevels(gb))
    seqnames(gb)[i] <- "chr0"
  } else {
    seqlevels(gb) <- c(seqlevels(gb), "chrZ")
    seqnames(gb)[i] <- "chrZ"
  }
  gb
}
