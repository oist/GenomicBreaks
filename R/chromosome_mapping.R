#' Chromosome Mapping
#' 
#' This function computes a mapping between the chromosomes of the reference and query genomes. 
#' The mapping is chosen to maximize the number of aligned blocks in the matched chromosomes. 
#' The function tries to build a \code{1:1} mapping, but this outcome is not guaranteed (see note below).
#' 
#' @param gb A [`GBreaks`] object.
#' @param min_block_size_perc A value between 0 and 100 used to filter the aligned blocks based on their size. 
#'                            Only aligned blocks with sizes greater than the percentile specified by the parameter are considered for mapping.
#'                            For example, \code{min_block_size_perc=10} will ignore the 10% smallest blocks.
#' 
#' @return Returns a list where each element corresponds to a mapped pair of chromosomes.
#'         Each element contains the following:
#'         1. \code{ref_chr_name} : the name of the reference chromosome;
#'         2. \code{qry_chr_name} : the name of the query chromosome;
#'         3. \code{gb} : a _GenomicBreaks_ object containing only aligned blocks for the corresponding chromosome
#'                   pair that are larger than the threshold specified by \code{min_block_size_perc}.
#'
#' @note The function tries to build a \code{1:1} mapping between chromosomes in the reference and query genomes; 
#'       however, in exceptional cases, one chromosome may end up mapped to multiple chromosomes. 
#'       When the genomes have unequal chromosome counts, the genome with more chromosomes will necessarily have unpaired chromosomes.
#'
#' @author Priscila Biller
#'
#' @family modifier functions
#'
#' @examples
#' \dontrun{
#' 	chrom_map <- chromosome_mapping(gb, 10)
#' 	cat("Chromosome mapping:\n")
#' 	for (chrom_pair in chrom_map) {
#' 	   cat(paste("- " , chrom_pair$ref_chr_name, " <-> ", chrom_pair$qry_chr_name, " (Nb. aligned blocks: ", length(chrom_pair$gb@ranges@start), ") \n", sep=""))
#' 	}
#' 	}
#'
#' @export
chromosome_mapping <- function(gb, min_block_size_perc=0){

	ref_chr_lst <- sort(levels(gb@seqnames@values))
	qry_chr_lst <- sort(levels(gb@elementMetadata@listData$query@seqnames@values))

	chr_lst_rows <- ref_chr_lst  # Rows: It refers to the genome with the least nb. of chromosomes.
	chr_lst_cols <- qry_chr_lst  # Columns: It refers to the genome with the most nb. of chromosomes.
	if(length(ref_chr_lst) > length(qry_chr_lst)){
		# Flipping reference and query genomes in the matrix: 
		# Reference chromosomes are now columns, and query chromosomes are now rows, 
		# as the query genome has less chromosomes than the reference.
		chr_lst_rows <- qry_chr_lst
		chr_lst_cols <- ref_chr_lst
	}

	# Create a map between chromosomes.
	M_chr <- matrix(0, nrow = length(chr_lst_rows), ncol = length(chr_lst_cols))
	for (i in seq_along(chr_lst_rows)) {
		row_chr <- chr_lst_rows[i]
		for (j in seq_along(chr_lst_cols)) {
			col_chr <- chr_lst_cols[j]
			# Get reference and query chromosomes.
			ref_chr <- row_chr
			qry_chr <- col_chr
			if(length(ref_chr_lst) > length(qry_chr_lst)){
				ref_chr <- col_chr
				qry_chr <- row_chr
			}
			# Compute number of aligned blocks.
			gb_filtered <- gb[seqnames(gb) == ref_chr & seqnames(gb$query) == qry_chr, ] |> coalesce_contigs()
			percentile_value <- quantile(width(gb_filtered), probs=min_block_size_perc/100)
			gb_filtered <- gb_filtered[width(gb_filtered) > percentile_value]
			N  <- length(gb_filtered@ranges@start)
			M_chr[i, j] <- N
		}
	}

	# Chromosome mapping: Column indices with highest values in each row.
	max_indices <- apply(M_chr, 1, which.max)

	# Get pairs of mapped chromosomes.
	chromosome_pairs <- list()
	for (i in seq_along(chr_lst_rows)) {

		max_j   <- max_indices[i]

		ref_chr <- ""
		qry_chr <- ""
		if(length(ref_chr_lst) > length(qry_chr_lst)){
			ref_chr <- ref_chr_lst[max_j]
			qry_chr <- qry_chr_lst[i]
		} else {
			ref_chr <- ref_chr_lst[i]
			qry_chr <- qry_chr_lst[max_j]
		}

		gb_filtered <- gb[seqnames(gb) == ref_chr & seqnames(gb$query) == qry_chr, ] |> coalesce_contigs() #|> (\(x) x[width(x) > min_block_size])()
		percentile_value <- quantile(width(gb_filtered), probs=min_block_size_perc/100)
		gb_filtered <- gb_filtered[width(gb_filtered) > percentile_value]

		N  <- length(gb_filtered@ranges@start)

		# Save info in a list where each element corresponds to a mapped pair of chromosomes.
		# An element in the list is also a list, composed of 3 elements: 
		# (1) Name of the chromosome in the reference genome;
		# (2) Name of the chromosome in the query genome;
		# (3) A GBreaks object.
		chromosome_pairs[[i]] <- list(ref_chr_name=ref_chr, qry_chr_name=qry_chr, gb=gb_filtered)
	}
	chromosome_pairs
}
