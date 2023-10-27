#' Index measuring how 'random' the alignment strand has become.
#'
#' In groups of species where major changes of gene order happened but genes
#' tended to stay on the same chromosome, a marked feature is that the strand
#' on which homologous sequences align appears to be randomised.  This index
#' expresses it with a number.
#'
#' For each feature of the _target_ genome, the total length of its alignments
#' on the plus strand is subtracted from the total length of its alignments
#' on the minus strand.  The sign of this number is removed because there is
#' no guarantee that homologous chromosomes are sequenced in the same direction
#' in different assemblies.  To make features comparable, this number is then
#' normalised by total number of aligned bases.  Then, the computed numbers for
#' all the _target_'s features are averaged with a weighted mean, where the
#' weight is the total number of aligned bases of that feature.  Thus, a number
#' close to 1 is expected for closely related genomes.
#'
#' @note This index is designed for comparison of chromosomal assemblies that
#' have a strong conservation of synteny in the sense most homologous genes
#' are mapped on pairs of homologous chromosomes.
#'
#' @param gb A [`GBreaks`] object.
#'
#' @family Similarity indexes
#'
#' @returns Returns a numeric value between 0 and 1.
#'
#' @examples
#' strand_randomisation_index(exampleInversion)
#'
#' @export

strand_randomisation_index <- function(gb) {
  if(length(gb) == 0) return(numeric(0))
  if(length(gb) == 1) return(    1     )
  gbl <- split(gb, seqnames(gb), drop = TRUE)
  # Calculate an index for each sequence feature
  idx <- sapply(gbl, \(x) {
    onPlus  <- sum(width(x[strand(x) == '+']))
    onMinus <- sum(width(x[strand(x) == '-']))
    abs((onPlus - onMinus) / (onPlus + onMinus))
  })
  # Average by the sum of all widths
  weighted.mean(idx, sum(width(gbl)))
}
