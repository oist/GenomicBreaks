#' Index measuring how 'random' the alignment strand has become.
#'
#' In groups of species where major changes of gene order happened but genes
#' tended to stay on the same chromosome, a marked feature is that the strand
#' on which homologous sequences align appears to be randomised.  This index
#' expresses it with a number.
#'
#' For each feature in the _target_ genome, the total length of alignments on
#' the plus strand is subtracted from that on the minus strand, and the absolute
#' value is taken (because there is no guarantee that homologous chromosomes are
#' sequenced in the same direction in different assemblies).  This value is
#' normalized by the total number of aligned bases.  A weighted mean is then
#' computed across all features, using each feature’s total aligned bases as its
#' weight.  Thus, a number close to 1 is expected for closely related genomes.
#'
#' @note This index is designed for comparison of chromosomal assemblies that
#' have a strong conservation of synteny in the sense most homologous genes
#' are mapped on pairs of homologous chromosomes.  In 2025, a new `tiles` option
#' was added to make the index more robust to assembly errors and karyotype
#' changes.
#'
#' @param gb A [`GBreaks`] object.
#' @param tiles A number of tiles
#'
#' @family Similarity indexes
#'
#' @returns Returns a numeric value between 0 and 1.
#'
#' @references Plessy C, Mansfield MJ, Bliznina A, Masunaga A, West C, Tan Y,
#' Liu AW, Grašič J, Del Río Pisula MS, Sánchez-Serna G, Fabrega-Torrus M,
#' Ferrández-Roldán A, Roncalli V, Navratilova P, Thompson EM, Onuma T, Nishida
#' H, Cañestro C, Luscombe NM (2024). Extreme genome scrambling in marine
#' planktonic _Oikopleura dioica_ cryptic species. *Genome Research*, 34(3),
#' 426–440. \url{https://doi.org/10.1101/gr.278295.123}. PMID: 38621828
#'
#' @author Charles Plessy
#'
#' @examples
#' strand_randomisation_index(exampleInversion)
#'
#' @export

strand_randomisation_index <- function(gb, tiles=NULL) {
  if(length(gb) == 0) return(numeric(0))
  if(length(gb) == 1) return(    1     )
  # Use the original sequence features if `tiles` is `NULL`
  if (is.null(tiles)) {
    gbl <- split(gb, seqnames(gb), drop = TRUE)
  } else {
    gb <- forceSeqLengths(gb)
    # Else, tile the genome
    tiles <- tile( GRanges(seqinfo(gb)), n=tiles) |> unlist() |> unname()
    # Create an object representing boundaries between tiles on both strands
    boundaries  <- narrow(tiles, 2) |> gaps() |> rep(2) |> sort()
    strand(boundaries) <- c("+", "-")
    # Cut the target genome sequences at each boundary
    cutByTiles <- setdiff(gb, boundaries) |> sort(ignore.strand = TRUE)
    gbl <- split(cutByTiles, subjectHits(findOverlaps(cutByTiles, tiles)))
  }
  # Calculate an index for each sequence feature
  idx <- sapply(gbl, \(x) {
    onPlus  <- sum(width(x[strand(x) == '+']))
    onMinus <- sum(width(x[strand(x) == '-']))
    abs((onPlus - onMinus) / (onPlus + onMinus))
  })
  # Average by the sum of all widths
  weighted.mean(idx, sum(width(gbl)))
}
