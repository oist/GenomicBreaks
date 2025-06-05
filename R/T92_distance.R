#' Index representing the evolutionary distance between two genomes.
#'
#' The T92 distance is calculated as -h * log(1 - p / h - q) - 0.5 * (1 - h) * log(1 - 2 * q). See [WIKI](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#T92_model_(Tamura_1992)) for more details.
#' @param list A list containing the probabilities of the alignment.
#'
#' @note Here, the probability ranges from 0 to 1. The list input can contain other information, but the function will only use the name in this format probability_{REF}_{ALT}
#' 
#' @family Similarity indexes
#'
#' @returns Returns a numeric value show the evolutionary distance between two genomes. the larger the value, the more different the two genomes are.
#'
#' @examples
#' parameters <- readTrainFile("train.txt")
#' T92_distance(parameters)
#'
#' @export

# T92 model (Tamura 1992)
T92_distance <- function(list) {

  theta <- prob_matrix[["G_G"]] + prob_matrix[["C_C"]]
  h <- 2 * theta * (1 - theta)
  p <- prob_matrix[["A_G"]] + prob_matrix[["G_A"]] + prob_matrix[["T_C"]] + prob_matrix[["C_T"]]
  q <- prob_matrix[["A_C"]] + prob_matrix[["C_A"]] + prob_matrix[["G_T"]] + prob_matrix[["T_G"]] + prob_matrix[["A_T"]] + prob_matrix[["T_A"]] + prob_matrix[["C_G"]] + prob_matrix[["G_C"]]
  d <- - h * log(1 - p / h - q) - 0.5 * (1 - h) * log(1 - 2 * q)
  return(d)

}
