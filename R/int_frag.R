#' Integrated intersection of both query and reference subsets on a continuous alignment
#'
#' @param og_r GRanges object containing the refernce genome information of the original alignment
#' @param og_q GRanges object containing the query genome information of the original alignment
#' @param sub_r GRanges object containing the reference subset
#' @param sub_q GRanges object containing the query subset
#' @return GRanges object of alignment intersected on both reference and query side. Query data held in metadata column
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb

int_frag <- function(og_r, og_q, sub_r, sub_q){
  if(isEmpty(intersect(og_r, sub_r)) | isEmpty(intersect(og_q, sub_q))){
    return(NULL)
  }
  else{
    # coordinates (changing ref -> query (ref - co_shift) and query -> ref (query + co_shift))
    co_shift <- start(og_r) - start(og_q)

    # ref intersection
    int_r <- intersect(og_r, sub_r)

    # query intersection
    int_q <- intersect(og_q, sub_q)

    # convert query intersection points to ref coords (construct gr object in ref coords)
    converted_q_start <- start(int_q) + co_shift
    converted_q_end <- end(int_q) + co_shift
    q_frags <- length(int_q) #number of intersected query fragments
    converted_q <- GRanges(seqnames = rep(seqnames(og_r), q_frags), IRanges(start = converted_q_start, end = converted_q_end))

    # intersect intersected ref and converted intersected query
    fin_r <- intersect(int_r, converted_q)
    r_frags <- length(fin_r)

    # covert fin_r to query coords
    fin_q_start <- start(fin_r) - co_shift
    fin_q_end <- end(fin_r) - co_shift
    fin_q <- GRanges(seqnames = rep(seqnames(og_q), r_frags), IRanges(start = fin_q_start, end = fin_q_end))

    # return final gr object
    fin_r$name <- fin_q
    return(fin_r)
  }


}
