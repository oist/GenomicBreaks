#' Index representing how the karyotype changes between chromosomes of two genomes
#'
#' The index is calculated as the mean of the number of chromosomes in the two genomes.  ###############
#' 
#'
#' @param gb A [`GBreaks`] object.
#'
#' @note Here, the definition of \sQuote{karyotype change} is 1 - (the mean value of change score for each chromosome in the reference genome. 
#' The change score for chromosome Ref is calculated as 1 / (number of mapped chromosomes in the query genome). To be noticed, 
#' the mapped length of mapped chromosome in the query genome should be at least 10% of the length of the reference chromosome.
#' if chromosome Ref has only one mapped chromosome in the query genome, the change score is 1. For 2 chromosomes, the change score is 0.5.
#' For the symmetry, the sum of change scores for each chromosome in query will be normalized to 1 if the sum is more than 1.
#' 
#' @family Similarity indexes
#'
#' @returns Returns a numeric value between 0 and 1. 0 is expected for same karyotype genomes.
#'
#' @examples
#' gb       <- GRanges(c("Ref:100-200:+",   "Ref:400-500:+",    "Re2:600-700:+"))
#' gb$query <- GRanges(c("Que:1100-1200:+", "Que2:1700-1800:+", "Que2:1500-1600:+"))
#' karyotype_change_index(gb)
#'
#' @export

karyotype_change_index <- function(gb) {

  if(length(gb) == 0) return(numeric(1))
  if(length(gb) == 1) return(    0     )
  gbl <- split(gb, seqnames(gb), drop = TRUE)
  
  # Calculate an index for each sequence feature
  chrScore_list <- lapply(gbl, function(x) {
    # Keep only matches to the best query
    allQuery <- tapply(width(x$query), seqnames(x$query), function(x) sum(x, na.rm = TRUE))
    
    # # filter by length of mapped sequence, can not deal with chrY
    # allMappedLength <- sum(allQuery, na.rm = TRUE)
    # mainQuery <- allQuery[!is.na(allQuery) & allQuery >= allMappedLength * 0.1]
    
    # filter by length of estimate chromosome length
    range_df <- as.data.frame(ranges(x))
    chr_length <- max(range_df$end) - min(range_df$start)
    mainQuery <- allQuery[!is.na(allQuery) & allQuery >= chr_length * 0.1]
    
    # final scores
    if (length(mainQuery) == 0) return(NULL) # return NULL or 0? need to explore
    ###print(mainQuery)
    ###print(length(mainQuery))
    ###print(rep(as.character(seqnames(x)[1]), length(mainQuery)))
    ###print(names(mainQuery))
    ###print(rep(1 / length(mainQuery), length(mainQuery)))
    data.frame(target = rep(as.character(seqnames(x)[1]), length(mainQuery)),
               query = names(mainQuery),
               value = rep(1 / length(mainQuery), length(mainQuery)), 
               stringsAsFactors = FALSE
                )
  })

  # filter out NULL
  chrScore_list <- chrScore_list[!sapply(chrScore_list, is.null)]

  # combine all data frames into one
  chrScore_df <- do.call(rbind, chrScore_list)
  rownames(chrScore_df) <- NULL
  
  chrScore_df <- chrScore_df %>%
    group_by(query) %>%
    mutate(value = case_when(
      sum(value) > 1 ~ value / sum(value, na.rm = TRUE),
      TRUE ~ value
    ))
  ###print(chrScore_df)
  
  #return(chrScore_df$value)
  return(1 - round(mean(chrScore_df$value), 4))

}
