#' Index representing how the karyotype changes between chromosomes of two genomes
#'
#' This function is used to read the train file from LAST and extract the final parameters.
#'
#'
#' @param input_file A string of the path to the train file from LAST.
#'
#'
#' @family Similarity indexes
#'
#' @returns Returns a list of parameters from LAST train file. "PercentSimilarity", "PercentSimilarityNogaps", "mean_delete_size", "mean_insert_size", "substitution_percent_identity", "probability",
#' "matchProb", "delOpenProb", "insOpenProb", "delExtendProb", "insExtendProb", "endProb"
#'
#' @examples
#' parameters <- readTrainFile(system.file("extdata/example.train", package = "GenomicBreaks"))
#'
#' @export

readTrainFile <- function(input_file) {

  extract_parameters <- function(record) {
    parts <- strsplit(record, ": ")[[1]]
    par <- list()
    key <- gsub(" ", "_", parts[1])
    values <- strsplit(parts[2], " ")[[1]]
    values <- values[values != ""]

    if (length(values) == 1) {
      par[[key]] <- as.numeric(values[1])
    } else {
      bases <- c("AT", "CG")
      for (i in seq_along(bases)) {
        subkey <- gsub("%", "", key)
        par[[paste0(subkey, bases[i])]] <- as.numeric(values[i])
      }
    }
    return(par)
  }

  # Function to extract matrix values
  extract_matrix <- function(records, prefix) {
    matrix <- list()
    header <- strsplit(records[1], "\\s+")[[1]]
    print(header)
    for (i in 2:5) {
      row <- strsplit(records[i], "\\s+")[[1]]
      row_name <- row[1]
      values <- as.numeric(row[2:length(row)])
      for (j in seq_along(values)) {
        key <- paste0(prefix, "_", row_name, "_", header[j + 1])
        matrix[[key]] <- values[j]
      }
    }
    return(matrix)
  }

# Read and process file
if (! file.exists(input_file)) {
  stop("File does not exist")
}

lines <- readLines(input_file)
records <- character()
collect <- FALSE

for (line in lines) {
  if (grepl("^# lastal", line)) {
    records <- character()
    collect <- TRUE
    next
  }
  if (startsWith(line, "#") && collect) {
    cleaned <- gsub("^# ", "", line)
    records <- c(records, cleaned)
  }
}

  # Parse records
  parameters <- list()
  i <- 1
  while (i <= length(records)) {
    record <- records[i]
    if (grepl("matrix", record)) {
      prefix <- strsplit(record, "\\s+")[[1]][1]
      matrix <- extract_matrix(records[(i+1):(i+5)], prefix)
      if (grepl("probability", names(matrix)[1])) {
        prob_matrix <- matrix
      }
      parameters <- c(parameters, matrix)
      i <- i + 5
    } else if (grepl(":", record)) {
      param <- extract_parameters(record)
      for (key in names(param)) {
        if (! key %in% names(parameters)) {
          parameters <- c(parameters, param[key])
        }
      }
    } else if (grepl("#last", record)) {
      break
    }
    i <- i + 1
  }

  keywords <- c("PercentSimilarity", "PercentSimilarityNogaps", "mean_delete_size", "mean_insert_size", "substitution_percent_identity", "probability",
                "matchProb", "delOpenProb", "insOpenProb", "delExtendProb", "insExtendProb", "endProb")
  filtered_parameters <- parameters[grepl(paste(keywords, collapse = "|"), names(parameters))]

}
