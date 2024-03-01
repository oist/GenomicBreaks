#install.packages("zoo")
#library(zoo)

pattern_search <- function(seq){

  diff_seq <- diff(seq)

  simple <- c(2, -1, 2)
  double <- c(5, 1, -4, -1, 3, -1, 4)
  double2 <- c(4, -1, 3, -1, -4, 1, 5)
  nested <- c(6, -1, -2, 1, -2, -1, 6)
  nested2 <- c(3,  1, -2, -1,  4)
  nested3 <- c(4, -1, -2,  1,  3)
  sequential <- c(2, -1, 3, -1, 2)
  sequential2 <- c(2, -1, 3, -1, 3, -1, 2)

  inversions <- matrix(NA, nrow=0, ncol=2)

  pattern <- simple
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+1, i+2))
  }

  pattern <- simple*(-1)
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+1, i+2))
  }

  pattern <- double
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+1, i+4)) |> rbind(c(i+3, i+6))
  }

  pattern <- double2
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+3, i+6)) |> rbind(c(i+1, i+4))
  }

  pattern <- double*(-1)
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+1, i+4)) |> rbind(c(i+3, i+6)) #left first, then right
  }

  pattern <- double2*(-1)
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+3, i+6)) |> rbind(c(i+1, i+4)) #right first, then left
  }

  pattern <- nested
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+1, i+6)) #simple will be find above
  }

  pattern <- nested*(-1)
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+1, i+6)) #simple will be find above
  }

  pattern <- nested2
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+1, i+2)) |> rbind(c(i+1, i+4)) #left first, then all
  }

  pattern <- nested2*(-1)
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+1, i+2)) |> rbind(c(i+1, i+4)) #left first, then all
  }

  pattern <- nested3
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+1, i+4)) |> rbind(c(i+1, i+2)) #all first, then left
  }

  pattern <- nested3*(-1)
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+1, i+4)) |> rbind(c(i+1, i+2)) #all first, then left
  }

  pattern <- sequential
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+1, i+2)) |> rbind(c(i+3, i+4))
  }

  pattern <- sequential*(-1)
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+1, i+2)) |> rbind(c(i+3, i+4))
  }

  pattern <- sequential2
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+1, i+2)) |> rbind(c(i+3, i+4)) |> rbind(c(i+5, i+6))
  }

  pattern <- sequential2*(-1)
  pattern_length <- length(pattern)
  matches <- which(rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left"))
  for (i in matches){
    inversions <- rbind(inversions, c(i+1, i+2)) |> rbind(c(i+3, i+4)) |> rbind(c(i+5, i+6))
  }

return(inversions)
}

find_unique_sequences <- function(numbers) {

  #generate k-length sequences and their frequencies

  numbers <- diff(numbers)

  generate_and_count_sequences <- function(numbers, k) {
    numbers_as_str <- sapply(numbers, as.character)
    sequences <- vector("list", length(numbers) - k + 1)
    for (i in 1:(length(numbers) - k + 1)) {
      sequences[[i]] <- paste("(", paste(numbers_as_str[i:(i+k-1)], collapse = " "), ")", sep="")
    }
    frequencies <- table(unlist(sequences))
    return(frequencies)
  }

  k <- 3
  all_patterns <- list()

  repeat {
    frequencies <- generate_and_count_sequences(numbers, k)
    all_patterns[[as.character(k)]] <- sort(frequencies, decreasing = TRUE)

    if (all(frequencies == 1)) {
      break
    }

    k <- k + 1
  }

  patterns_different_from_one <- lapply(all_patterns, function(f) {
    f[f > 5]
  })

  return(patterns_different_from_one)
}
