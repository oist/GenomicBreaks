#' Variation of the Jukes-Cantor 1969 distance
#'
#' The Jukes-Cantor 1969 (JC69) distance corrects the p-distance for multiple substitutions, 
#' providing an estimate of evolutionary distance that is proportional to time under the model.
#' The correction is based on the proportion of nucleotide differences, typically obtained 
#' by counting mismatches in the alignment matrix.
#' Here, the Jukes-Cantor equation remains unchanged, but the definition of which base pairs 
#' are considered different is modified.
#'
#' In this function, the fraction of nucleotides that are different incorporates not only the 
#' mismatches from the alignment matrix but also the base pairs that are left unaligned.
#' The rationale is that these unaligned base pairs likely differ primarily because of point 
#' substitutions and should therefore be treated as mismatches that were not detected by the aligner. 
#' In fact, the aligner can only spot mismatches in regions where the proportion of matches is high enough for alignment.
#'
#' Notice that gaps are usually not included in the Jukes-Cantor distance, including in this variation, 
#' because they are generally considered to result from indels (insertions and deletions that affect multiple 
#' nucleotides at once), whereas the Jukes-Cantor model is based only on point substitutions. 
#' Therefore, including regions that were likely affected by large evolutionary events, such as gaps caused by indels, 
#' would incorrectly inflate a distance calculated under a model in which only one position is mutated at a time.
#'
#' @references Jukes, T.H. & Cantor, C.R. (1969). "Evolution of protein molecules." In *Mammalian Protein Metabolism* (pp. 21–132). Academic Press.
#'
#' @param gb A [`GBreaks`] object.
#' @param m  A matrix of **counts** for bases of the _target_ genome to be aligned to bases on the _query_ genome.
#' @param adjust_p A boolean flag. If `TRUE`, the distance is scaled between `0` and `0.75` to ensure the logarithm stays positive.
#'
#' @family Alignment statistics
#' @family Similarity indexes
#'
#' @author Priscila Biller
#'
#' @returns Returns a numeric value representing the evolutionary distance between two genomes. The greater the value, the more genetically different the genomes are.
#'
#' @examples
#'
#' # Only the sequence length is used from the GenomicBreaks object.
#' gb       <- GRanges(c("Ref:100-35000000:+"))
#' gb$query <- GRanges(c("Que:1100-35000500:+"))
#' d <- JC69_distance_allseq(gb, exampleSubstitutionMatrix)
#'
#' @export
JC69_distance_allseq <- function(gb, m, adjust_p=FALSE) {
  
  if(length(gb) == 0) return(numeric(0))
  if (all(m == 0)) return(NA)

  # Gets the smallest sequence length.
  totBps <- min(sum(guessSeqLengths(gb)),sum(guessSeqLengths(gb$query)))

  # Matrix of aligned base pairs, excluding gaps.
  non_gap    <- c("A", "C", "G", "T")
  m_non_gap  <- m[non_gap, non_gap]

  aligned    <- sum(m)                     # Total aligned base pairs.
  matches    <- sum(diag(m_non_gap))       # Matches.
  mismatches <- sum(m_non_gap)-matches     # Mismatches.
  gaps       <- aligned-matches-mismatches # Gaps
  unaligned  <- totBps - aligned           # Total unaligned base pairs.

  JC69_distance(mismatches+unaligned, tot=totBps-gaps, adjust_p=adjust_p)
} 
