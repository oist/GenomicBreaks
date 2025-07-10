#' Jukes-Cantor 1969 distance
#'
#' The Jukes and Cantor 1969 (JC69) distance is calculated as \eqn{-\frac{3}{4} \ln(1 - \frac{4}{3} p )}, where \eqn{p} is the sum of the probabilities of the alignment. See the [Wikipedia](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_(Jukes_and_Cantor_1969)) for more details.
#'
#' @references Jukes, T.H. & Cantor, C.R. (1969). "Evolution of protein molecules." In *Mammalian Protein Metabolism* (pp. 21–132). Academic Press.
#'
#' @param train_parameters A list containing the probabilities of the alignment, produced by the [`readTrainFile()`] function.
#'
#' @family Similarity indexes
#'
#' @author Zikun Yang
#'
#' @returns Returns a numeric value show the evolutionary distance between two genomes. the larger the value, the more different the two genomes are.
#'
#' @examples
#' parameters <- readTrainFile(system.file("extdata/example.train", package = "GenomicBreaks"))
#' JC69_distance(parameters)
#'
#' @export

JC69_distance <- function(train_parameters) {
  p <- train_parameters[["probability_A_C"]] + train_parameters[["probability_A_G"]] + train_parameters[["probability_A_T"]] + train_parameters[["probability_C_A"]] + train_parameters[["probability_C_G"]] + train_parameters[["probability_C_T"]] + train_parameters[["probability_G_A"]] + train_parameters[["probability_G_C"]] + train_parameters[["probability_G_T"]] + train_parameters[["probability_T_A"]] + train_parameters[["probability_T_C"]] + train_parameters[["probability_T_G"]]
  - 0.75 * log(1 - 4 * p / 3)
}
