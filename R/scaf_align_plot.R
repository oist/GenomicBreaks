#' Plot Alignment of Single Reference Scaffold
#'
#' The function plots the pairwise alignments from one of the reference scaffolds
#'
#' @param gr_ob GRanges object containing pairwise alignment
#' @param scaf Character string of the name of the scaffold you wish to plot
#' @return A plot of the query scaffold alignments to the chosen reference scaffold
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom graphics plot.new plot.window rect


scaf_align_plot <- function(gr_ob, scaf){

  # Prep
  ref_now <- gr_ob[seqnames(gr_ob) == scaf]
  q_now <- ref_now$query

  q_scafs <- as.numeric(seqnames(q_now))

  if (all(is.na(seqlengths(gr_ob)))){
    ref_len <- max(end(gr_ob))
  }
  else{
    ref_len <- seqlengths(gr_ob)[scaf]
  }

  s_len <- length(ref_now)

  # Plotting
  plot.new()
  plot.window(xlim = c(0, ref_len), ylim = c(min(q_scafs), max(q_scafs)*1.1 +0.5))
  title(xlab = paste("Reference Genome Scaffold", scaf), ylab = "Scaffolds of Query Genome")

  axis(1)
  axis(2)

  rect(xleft = start(ranges(ref_now)), xright = end(ranges(ref_now)),
       ybottom = q_scafs - 0.05, ytop = q_scafs + 0.05)

  # Coverage
  rect(xleft = start(ranges(ref_now)), xright = end(ranges(ref_now)),
       ybottom = as.vector(rep(max(q_scafs)*1.1 -0.05, s_len)),
       ytop = as.vector(rep(max(q_scafs)*1.1 +0.05, s_len)))

}
