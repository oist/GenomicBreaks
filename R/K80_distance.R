#' Index representing the evolutionary distance between two genomes.
#'
#' The K80 distance is calculated as -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q)). See [WIKI](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#K80_model_(Kimura_1980)) for more details.
#' 
#' @param list A list containing the probabilities of the alignment.
#'
#' @note Here, the probability ranges from 0 to 1. The list input can contain other information, but the function will only use the name in this format probability_{REF}_{ALT}
#' 
#' @family Similarity indexes
#'
#' @returns Returns a numeric value show the evolutionary distance between two genomes. the larger the value, the more different the two genomes are.ss
#'
#' @examples
#' parameters <- readTrainFile("train.txt")
#' K80_distance(parameters)
#'
#' @export

# K80 model (Kimura 1980)
K80_distance <- function(list) {

  p <- prob_matrix[["A_G"]] + prob_matrix[["G_A"]] + prob_matrix[["T_C"]] + prob_matrix[["C_T"]]
  q <- prob_matrix[["A_C"]] + prob_matrix[["C_A"]] + prob_matrix[["G_T"]] + prob_matrix[["T_G"]] + prob_matrix[["A_T"]] + prob_matrix[["T_A"]] + prob_matrix[["C_G"]] + prob_matrix[["G_C"]]
  k <- - 0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q))
  return(k)

}
