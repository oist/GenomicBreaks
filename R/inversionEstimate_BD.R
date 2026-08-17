
#' Expected number of cycler after k reversals
#' 
#' Computes the expected number of cycles in a random graph
#' after k edges. 
#' 
#' Implements Theorem 3 from Berestycki and Durrett (2006),
#' also described in Equation 6 in Biller et al. (2015).
#'
#' @param k An integer: The number of edges in the random graph, 
#'                      equivalent to the number of inversions.
#'
#' @param N An integer: The number of vertices in the random graph, 
#'                      equivalent to the number of markers+1.
#' 
#' @references Berestycki, Nathanaël, and Rick Durrett. "A phase transition in the random transposition random walk." Discrete Mathematics and Theoretical Computer Science. Discrete Mathematics and Theoretical Computer Science, 2003.
#' @references Hannenhalli, Sridhar, and Pavel A. Pevzner. "Transforming cabbage into turnip: polynomial algorithm for sorting signed permutations by reversals." Journal of the ACM (JACM) 46.1 (1999): 1-27.
#' @references Biller, Priscila, Laurent Guéguen, and Eric Tannier. "Moments of genome evolution by double cut-and-join." BMC bioinformatics 16.Suppl 14 (2015): S7.
#'
#' @author Priscila Biller
#'
#' @export
#' @keywords internal
nb_cycles_k <- function(k, N){

    # Invalid value (negative number of inversions or not enough markers).
    if ((N < 2) | (k < 0)){
      NA

    # No inversions. Genomes are the same.
    } else if(k == 0) {
      N

    # At least one inversion.
    } else {

      nb_cyc_k <- 0
      # This loop is supposed to go until infinity (sum until infinity), 
      # but as i grows, the contribution of the i-th term gets smaller and smaller.
      # Here, the computation is capped at 10 for efficiency. 
      # It does not seem to affect much the accuracy.
      for(i in 1:10) {
        term_base <- 2.0*k/N
        term_k    <- term_base*exp(-term_base)
        nb_cyc_k  <- nb_cyc_k + ((1.0 / term_base) * ( i**(i-2) / factorial(i) )* term_k**i)
      }
      nb_cyc_k*N
    }
}

#' Inversion Estimate - setup for option 'many'
#' 
#' The setup step computes the expected number of cycles 
#' after k edges for various values of k. 
#'
#' The expected number of cycles can take any value between 1 
#' and n (the number of markers). For each of these values, 
#' the method computes at least one k.
#' 
#' Notice that this is the most time-consuming step 
#' in the estimate. Since it depends only on the number
#' n of markers, it is recommended that it be computed once 
#' for all permutations of n markers, saving time.
#'
#' @param n An integer: The number of markers. 
#'          A marker represents an aligned block, a gene, 
#'          or any other genomic region being rearranged.
#'
#' @author Priscila Biller
#'
#' @export
#' @keywords internal
inversionEstimate_BD_setup <- function(n){

  # The components of the breakpoint graph for random reversals 
  # on n-1 markers are related to the cycles of random 
  # transposition on n markers.
  N <- n+1

  cyc_all <- list(list(k=0,nb_cyc_k=N))
  k       <- 1
  pos_k   <- 2

  while(cyc_all[[length(cyc_all)]]$nb_cyc_k > 1) {

    # Computes the expected number of cycles in a random graph after k edges.
    cyc_all[[pos_k]] <- list(k=k,nb_cyc_k=nb_cycles_k(k, N))

    # As the number of events increases, the expected number of cycles 
    # requires more events to be reduced by one. To avoid computing many values of k
    # that yield approximately the same expected number of cycles, the method skips some ks: 
    # it jumps when the difference in the number of cycles between the most recently 
    # computed value and its predecessor is less than one.
    k_diff      <- cyc_all[[pos_k]]$k - cyc_all[[pos_k-1]]$k
    # The number of cycles decreases as the number of events increases.
    nb_cyc_diff <- as.integer(cyc_all[[pos_k-1]]$nb_cyc_k) - as.integer(cyc_all[[pos_k]]$nb_cyc_k)

    if((nb_cyc_diff < 2) || (k_diff < 2)){
      last_pos    <- length(cyc_all)
      k_diff      <- cyc_all[[last_pos]]$k - cyc_all[[last_pos-1]]$k
      nb_cyc_diff <- as.integer(cyc_all[[last_pos]]$nb_cyc_k) - as.integer(cyc_all[[last_pos-1]]$nb_cyc_k)
      if(nb_cyc_diff == 0){
        k_diff <- k_diff*2
      }
      k     <- cyc_all[[last_pos]]$k + k_diff
      pos_k <- last_pos + 1
    } else {
      # Append a new value in the current position.
      k       <- cyc_all[[pos_k-1]]$k + as.integer(k_diff)/2
      cyc_all <- c(cyc_all[1:(pos_k-1)], NA, cyc_all[pos_k:length(cyc_all)])
    }
  }
  cyc_all
}

#' Inversion Estimate - option 'single'
#' 
#' It finds the number of inversions whose expected number of cycles 
#' is closest to the observed number of cycles.
#' 
#' This function is an alternative to the function "inversionEstimate_BD_many",
#' which uses the function "inversionEstimate_BD_setup".
#'
#' This function should be preferred when there are few inversion estimates to be computed.
#' 
#' @param n An integer: The number of markers. 
#'          A marker represents an aligned block, a gene, 
#'          or any other genomic region being rearranged.
#' 
#' @param obs_nb_cycles An integer: The observed number of cycles
#'                      in the breakpoint graph of two genomes.
#'
#' @author Priscila Biller
#'
#' @export
#' @keywords internal
inversionEstimate_BD_single <- function(n, obs_nb_cycles){

  # The components of the breakpoint graph for random reversals 
  # on n-1 markers are related to the cycles of random 
  # transposition on n markers.
  N <- n+1

  # At the start of the evolutionary process, when the two genomes are identical, 
  # the number of cycles is equal to the number of markers, the maximum possible value.
  # Knowing that an inversion can split one cycle into two or merge two cycles into one,
  # if c cycles are observed, then at least N-c inversions must have occurred.
  k_beg <- 0
  k_end <- N

  # These variables hold the previous k_beg and k_end (to avoid infinite loops):
  # If this case occurs, there is a bug.
  k_beg_prev <- k_beg-1
  k_end_prev <- k_end+1

  # These variables hold the last valid k_beg and k_end:
  # ------[k_beg ------ point ------ k_end]------
  k_beg_valid <- 0
  k_end_valid <- 100000*k_end

  nb_cyc_k_beg <- N
  nb_cyc_k_end <- nb_cycles_k(k_end, N)

  # Identical permutations. No inversions occurred.
  if(obs_nb_cycles == nb_cyc_k_beg){
    return(0)
  }

  #cat(paste("START: [", k_beg, " (", nb_cyc_k_beg, "), ", k_end, " (", nb_cyc_k_end, ")] / ", obs_nb_cycles, "\n"))
  while(!((as.integer(nb_cyc_k_beg) == obs_nb_cycles) & (as.integer(nb_cyc_k_end) == (obs_nb_cycles - 1))) & ((k_beg != k_beg_prev) | (k_end != k_end_prev))){

    # Update previous values.
    k_beg_prev <- k_beg
    k_end_prev <- k_end

    # Case: k_end does not encompass point
    # ------[k_beg ------ k_end]------ point ------
    if(as.integer(nb_cyc_k_end) >= obs_nb_cycles){

      k_beg_valid <- k_end

      k_beg <- k_end
      k_end <- k_end*2

      nb_cyc_k_beg <- nb_cyc_k_end
      nb_cyc_k_end <- nb_cycles_k(k_end, N)

    # Case: k_beg does not encompass point
    # ------point ------ [k_beg ------ k_end]------ 
    } else if(as.integer(nb_cyc_k_beg) < obs_nb_cycles){
      
      k_end_valid <- k_beg

      k_beg <- k_beg_valid + floor((k_beg-k_beg_valid)/2)
      k_end <- k_end_valid

      nb_cyc_k_end <- nb_cyc_k_beg
      nb_cyc_k_beg <- nb_cycles_k(k_beg, N)

    # Case: point is encompassed but interval is not tight on the left
    # ------[k_beg ------ point ------ k_end]------
    } else if(as.integer(nb_cyc_k_beg) > obs_nb_cycles){

      k_beg_valid <- max(k_beg, k_beg_valid)

      k_beg <- k_beg + floor((k_end-k_beg)/2)
      nb_cyc_k_beg <- nb_cycles_k(k_beg, N)

    # Case: point is encompassed but interval is not tight on the right
    # ------[k_beg point ------ k_end]------
    } else if(as.integer(nb_cyc_k_end) < (obs_nb_cycles-1)){
      k_end_valid <- min(k_end, k_end_valid)

      k_end <- k_beg + ceiling((k_end-k_beg)/2)
      nb_cyc_k_end <- nb_cycles_k(k_end, N)
    }

    #cat(paste(" [", k_beg, " (", nb_cyc_k_beg, "), ", k_end, " (", nb_cyc_k_end, ")] / ", obs_nb_cycles, "\n"))
  }

  #cat(paste("RES: [", k_beg, " (", nb_cyc_k_beg, "), ", k_end, " (", nb_cyc_k_end, ")] / ", obs_nb_cycles, "\n"))
  list(k_beg=k_beg, k_end=k_end, k_avg=as.integer((k_beg + k_end)/2), nb_cyc_beg=nb_cyc_k_beg, nb_cyc_end=nb_cyc_k_end)
}

#' Inversion Estimate - option 'many'
#' 
#' This function should be preferred when there are many inversion estimates to be computed.
#' 
#' @param n An integer: The number of markers. 
#'          A marker represents an aligned block, a gene, 
#'          or any other genomic region being rearranged.
#' 
#' @param obs_nb_cycles An integer: The observed number of cycles
#'                      in the breakpoint graph of two genomes.
#'
#' @author Priscila Biller
#'
#' @export
#' @keywords internal
inversionEstimate_BD_many <- function(n, obs_nb_cycles, cyc_all=NA){
  # Computes for various values of k the relation:
  # Number of inversions k <-> Expected number of cycles after k inversions
  if(is.na(cyc_all)){
    cyc_all <- inversionEstimate_BD_setup(n)
  }
  # Estimates expected number of inversions.
  all_nb_cycles <- sapply(cyc_all, '[[', "nb_cyc_k")
  closest  <- cyc_all[[which.min(abs(all_nb_cycles - obs_nb_cycles))]]
  list(k_beg=closest$k, k_end=closest$k, k_avg=closest$k, nb_cyc_beg=closest$nb_cyc_k, nb_cyc_end=closest$nb_cyc_k)
}

#' Inversion Estimate - Berestycki and Durrett (2006)
#' 
#' Computes the expected number of inversions that explains the differences between the gene orders of two genomes.
#' 
#' To compute the expected number of inversions, the method needs to first compute the *breakpoint graph*,
#' a classical data structure introduced by Hannenhalli and Pevzner and used in many genome 
#' rearrangement problems.
#' 
#' Given the number of cycles in the breakpoint graph,
#' it returns the estimated number of evolutionary events (inversions).
#' 
#' This function implements Theorem 3 from Berestycki and Durrett (2006),
#' also described in Equation 6 in Biller et al. (2015).
#' 
#' @references Berestycki, Nathanaël, and Rick Durrett. "A phase transition in the random transposition random walk." Discrete Mathematics and Theoretical Computer Science. Discrete Mathematics and Theoretical Computer Science, 2003.
#' @references Hannenhalli, Sridhar, and Pavel A. Pevzner. "Transforming cabbage into turnip: polynomial algorithm for sorting signed permutations by reversals." Journal of the ACM (JACM) 46.1 (1999): 1-27.
#' @references Biller, Priscila, Laurent Guéguen, and Eric Tannier. "Moments of genome evolution by double cut-and-join." BMC bioinformatics 16.Suppl 14 (2015): S7.
#'
#' @note If the input genomes are multichromosomal, it is recommended to call
#' the function \code{\link{matchPairs}} before this function, to create a mapping
#' between reference and query chromosomes.
#'
#' @param gb A [`GBreaks`] object.
#' @param mode Optional parameter that specifies how estimates are computed. Two options are available: 
#'             1. `single` : recommended when computing a few values in one function call; 
#'             2. `many` : use when computing many values. 
#' @param chrom_sep An optional string that separates the chromosome names in each item of the returned list.
#'
#' @return A list where each item corresponds to a chromosome pair. 
#'         The name of each item is the matched reference and query chromosome names concatenated, separated by `chrom_sep`.
#'         Each item contains the following information:
#'         1. Basic stats of the breakpoint graph: `N`, `nbBreakpoints`, `nbCycles`. See the documentation of `breakpointGraphProperties` for more details;
#'         2. Details of the _Berestycki and Durrett_ estimate stored in `expinv_BD`: 
#'            * `k_beg` and `k_end`: the lower and upper bounds for the expected number of inversions, respectively;
#'            * `k_avg` : the average of `k_beg` and `k_end`;
#'            * `nb_cyc_beg` and `nb_cyc_end`: the expected number of cycles in the breakpoint graph after `k_beg` and `k_end` inversions.
#'               The observed number of cycles in the breakpoint graph of the two genomes should fall between those values.
#'
#' @examples
#' \dontrun{
#' 
#' # Create a chromosome mapping given a GBreaks object (useful when genomes are multichromosomal).
#' chrMapping <- matchPairs(exampleInversionBader2001)
#' # Compute the expected number of inversions using the method from Berestycki and Durrett (2006).
#' expNbInversions <- inversionEstimate_BD(chrMapping)
#' # Compute the minimum number of inversions using the method from Hannehalli and Pevzner (1999).
#' minNbInversions <- inversionDistance(chrMapping)
#' # Output: Example from Bader et al. (2001): Inversion distance = 7 , Expected nb. of inversions = 7
#' cat(paste("Example from Bader et al. (2001): Inversion distance =", minNbInversions, ", Expected nb. of inversions =", expNbInversions[[1]]$expinv_BD$k_avg))
#' 
#' # Another example, this time without computing the chromosome mapping.
#' # The chromosome mapping is not needed if genomes are unichromosomal.
#' expNbInversions <- inversionEstimate_BD(exampleInversionBergeron2005b)
#' minNbInversions <- inversionDistance(exampleInversionBergeron2005b)
#' # Output: Example used in the book ``Mathematics of Evolution and Phylogeny`` (2005) (Figure 10.6): Inversion distance = 13 , Expected nb. of inversions = 15
#' cat(paste("Example used in the book ``Mathematics of Evolution and Phylogeny`` (2005) (Figure 10.6): Inversion distance =", minNbInversions, ", Expected nb. of inversions =", expNbInversions[[1]]$expinv_BD$k_avg))
#' }
#'
#' @seealso \code{\link{breakpointGraphProperties}} for computing key properties of the breakpoint graph needed for this estimate.
#'
#' @family Rearrangement distances
#' @family Similarity indexes
#'
#' @author Priscila Biller
#'
#' @export
inversionEstimate_BD <- function(gb, mode="single", chrom_sep="---"){

  # Get matched chromosomes.
  chrom_pairs <- lapply(split(gb, paste(seqnames(gb),seqnames(gb$query),sep=chrom_sep), drop=TRUE), function(x) {breakpointGraphProperties(x)})

  # Sort chromosomes from the one with the smallest number 
  # of aligned blocks to the one with the highest number.
  Ns <- sapply(chrom_pairs, `[[`, "N")
  chrom_pairs <- chrom_pairs[order(Ns)]

  N <- -1
  cyc_all <- NA

  for(chrom_pair in names(chrom_pairs)){
    chr <- chrom_pairs[[chrom_pair]]

    # If the chromosome has only one aligned block, there is nothing to be rearranged.
    if(chr$N == 1){
      chrom_pairs[[chrom_pair]]$expinv_BD <- NA
      next
    }

    # Mode: "many"
    if(mode == "many"){

      # If two chromosomes have the same number of markers,
      # then the most computationally expensive part of the 
      # estimate can be reused, saving time. 
      # Otherwise, a new computation is required.
      if(chr$N != N){
        cyc_all <- NA
        N <- chr$N
      }

      # Estimate expected number of inversions.
      chrom_pairs[[chrom_pair]]$expinv_BD <- inversionEstimate_BD_many(chr$N, chr$nb_cycles, cyc_all=cyc_all)

    # Mode: "single"
    } else {
      chrom_pairs[[chrom_pair]]$expinv_BD <- inversionEstimate_BD_single(chr$N, chr$nb_cycles)
    }
  }
  chrom_pairs
}
