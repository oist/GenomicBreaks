library(zoo)

pattern_search <- function(seq){

  diff_seq <- diff(seq)

  simple <- c(2, -1, 2)
  double <- c(5, 1, -4, -1, 3, -1, 4)
  double2 <- c(4, -1, 3, -1, -4, 1, 5)
  nested <- c(6, -1, -2, 1, -2, -1, 6)
  nested2 <- c(3,  1, -2, -1,  4)
  nested3 <- c(4, -1, -2,  1,  3)

  patterns <- list(simple, double, double2, nested, nested2, nested3)

  for (pattern in patterns){

    pattern_length <- length(pattern)
    matches <- rollapply(diff_seq, pattern_length, function(window) all(window == pattern), align = "left")
    print(which(matches))

    matches <- rollapply(diff_seq, pattern_length, function(window) all(window == pattern*(-1)), align = "left")
    print(which(matches))
  }

}
