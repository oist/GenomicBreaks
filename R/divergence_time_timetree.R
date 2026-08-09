#' Divergence time - TimeTree
#' 
#' Retrieves the divergence time of two species via the TimeTree API.
#' 
#' @param taxid1 An integer: The NCBI id of one taxon
#' @param taxid2 An integer: The NCBI id of the other taxon
#'
#' @return A list in which the first item is the median divergence time (\code{divtime_med_age}),
#'         and the last two items are the divergence time confidence interval 
#'         (\code{divtime_ci_low} and \code{divtime_ci_high}).
#'         If the divergence time is not found in the TimeTree database, it returns \code{NA}.
#'
#' @examples
#' \dontrun{
#' get_divergence_time(9555,9601)}
#'
#' @family Query functions
#'
#' @importFrom httr GET content
#'
#' @author Priscila Biller
#'
#' @export
get_divergence_time <- function(taxid1, taxid2) {
  timetree_pairwise_base <- "http://timetree.temple.edu/api/pairwise/"
  # Example of request to the TimeTree's API: 
  # - Request: http://timetree.temple.edu/api/pairwise/9555/9601
  # - Answer from TimeTree: 
  # taxon_a_id,taxon_b_id,scientific_name_a,scientific_name_b,all_total,precomputed_age,precomputed_ci_low,precomputed_ci_high,adjusted_age
  # 9555,9601,Papio anubis,Pongo abelii,83,28.82,26.8,30.6,0
  timetree_request <- paste(timetree_pairwise_base, taxid1, "/", taxid2, sep="")
  response <- GET(timetree_request)
  lines    <- content(response, as="text", encoding="UTF-8")
  # Break the whole content into lines, and gets the second line (data).
  line_data <- strsplit(lines, "\n", fixed = TRUE)[[1]][2] # First line: header; Second line: data.
  # Break the line into data.
  data <- strsplit(line_data, ",", fixed = TRUE)[[1]]
  # Check if the request returned something valid.
  if(length(data) >= 8) {
    list(divtime_med_age=as.numeric(data[6]),divtime_ci_low=as.numeric(data[7]),divtime_ci_high=as.numeric(data[8]))
  } else {
    NA
  }
}
