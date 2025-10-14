#' Tau index
#'
#' Tissue-specificity index from Yanai and coll., 2005 repurposed in 2022 by
#' Ranz and coll. for synteny analysis.
#'
#' In 2022, Ranz and coll., the tau index was computed for each chromosome in a
#' _target_ genome compared to a _query_, with a formula such as: for a given
#' _target_ feature, count its one-to-one orthologues on every feature \eqn{i}
#' of the _query_ genome, normalize these counts by dividing by the largest
#' count, subtract each normalized value from one, sum the results, and divide
#' by the number \eqn{n} of _query_ features minus one.
#'
#' \deqn{\frac{1}{n-1}\sum_{i=1}^{n}\left(1 - \frac{x_i}{\max(x)}\right)}{sum(1 - x / max(x)) / (length(x) - 1)}
#'
#' Here, the index is extended to the whole genome comparisons by computing it
#' for each feature of the _target_ genome and returning the average weighted
#' by feature length.
#'
#' @note Note that calculating the tau index on whole-genome nucleotide
#' alignments is not expected to produce meaningful results.  This function is
#' more useful when comparing the positions of protein orthologues.
#'
#' @param gb A [`GBreaks`] object.
#'
#' @family Similarity indexes
#'
#' @returns Returns a numeric value between 0 and 1.
#'
#' @references Yanai I, Benjamin H, Shmoish M, Chalifa-Caspi V, Shklar M, Ophir
#' R, Bar-Even A, Horn-Saban S, Safran M, Domany E, Lancet D, Shmueli O. (2005).
#' Genome-wide midrange transcription profiles reveal expression level
#' relationships in human tissue specification. *Bioinformatics*, 21(5):650-659.
#' \doi{10.1093/bioinformatics/bti042}. PMID: 15388519
#'
#' @references Ranz JM, González PM, Su RN, Bedford SJ, Abreu-Goodger C, Markow
#' T. (2022).  Multiscale analysis of the randomization limits of the
#' chromosomal gene organization between Lepidoptera and Diptera.
#' *Proc Biol Sci*, 289(1967):20212183. \doi{10.1098/rspb.2021.2183}.
#' PMID: 35042416
#'
#' @author Charles Plessy
#'
#' @examples
#' tau_index(exampleTranslocation)
#' GenomicBreaks:::.tau_index(c(0,8,0,0,0,2,0,2,0,0,0))
#'
#' @export

tau_index <- function(gb) {
  if(length(gb) == 0) return(numeric(0))
  if(length(gb) == 1) return(    1     )
  gbl <- split(gb, seqnames(gb), drop = TRUE)
  # Calculate an index for each sequence feature
  idx <- sapply(gbl, \(x) {
    gb$query |> seqnames() |> table() |> .tau_index()
  })
  # Average by the sum of all widths
  weighted.mean(idx, sum(width(gbl)))
}

.tau_index <- function(x) {
  if (length(x) == 0) return (numeric(0))
  if (length(x) == 1) return (1)
  sum ( 1 - x / max(x) ) / ( length (x) - 1)
}
