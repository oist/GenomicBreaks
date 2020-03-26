#' Compare Alignments over two Reference scaffolds
#'
#' The function plots the pairwise alignments from two of the reference scaffolds. This can be used as an aid in super-scaffolding. X-axis is the position on the reference scaffolds, and y-axis is aligned query scaffold.
#'
#' @param gr_ob GRanges object containing pairwise alignment
#' @param scaf A list containing two character objetcs of the names of the reference scaffolds you wish to compare on the plot
#' @return scaffold alignment plot, comparing two reference scaffolds
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom graphics plot.new plot.window rect axis title


align_scaf_plot_2<- function(gr_ob, scaf){


  ### reference scaffold 1 ###
  scaf1 <- scaf[[1]]
  ref_1 <- gr_ob[seqnames(gr_ob) == scaf1]
  q_1 <- GRanges(ref_1$name)

  q_scafs1 <- as.numeric(gsub("\\D+", "", as.vector(seqnames(q_1)))) #numeric value of query1 scaffolds

  ref1_len <- seqlengths(gr_ob)[scaf1]
  s_len1 <- length(end(ranges(ref_1)))

  ### reference scaffold 2 ###
  scaf2 <- scaf[[2]]
  ref_2 <- gr_ob[seqnames(gr_ob) == scaf2]
  q_2 <- GRanges(ref_2$name)

  q_scafs2 <- as.numeric(gsub("\\D+", "", as.vector(seqnames(q_2)))) #numeric value of query2 scaffolds

  ref2_len <- seqlengths(gr_ob)[scaf2]
  s_len2 <- length(end(ranges(ref_2)))

  ### combined info ###
  scaf_min <- min(min(q_scafs1), min(q_scafs2))
  scaf_max <- max(max(q_scafs1), max(q_scafs2))
  x_ax_max <- ref1_len + ref2_len

  ### Plotting ###
  plot.new()
  plot.window(xlim = c(0, x_ax_max), ylim = c(scaf_min, scaf_max*1.1 +0.05))
  title(xlab = paste("Ref Gen Scaffolds", gsub("\\D+", "", scaf[1]), "(left) and",
                     gsub("\\D+", "", scaf[2]), "(right)"), ylab = "Scaffolds of query Genome", main = "")
  axis(1)
  axis(2)

  ## Left Plot
  rect(xleft = start(ranges(ref_1)), xright = end(ranges(ref_1)),
       ybottom = q_scafs1 - 0.05, ytop = q_scafs1 + 0.05)

  ## Right Plot
  rect(xleft = start(ranges(ref_2)) + ref1_len, xright = end(ranges(ref_2)) + ref1_len, col = "blue", density = 100, border = "blue",
       ybottom = q_scafs2 - 0.05, ytop = q_scafs2 + 0.05)

  ## Coverage
  # Left
  rect(xleft = start(ranges(ref_1)), xright = end(ranges(ref_1)),
       ybottom = as.vector(rep(scaf_max*1.1 - 0.05, s_len1)),
       ytop = as.vector(rep(scaf_max*1.1 + 0.05, s_len1)))

  # Right
  rect(xleft = start(ranges(ref_2)) + ref1_len, xright = end(ranges(ref_2)) + ref1_len, col = "blue", density = 100, border = "blue",
       ybottom = as.vector(rep(scaf_max*1.1 - 0.05, s_len2)),
       ytop = as.vector(rep(scaf_max*1.1 + 0.05, s_len2)))
}

