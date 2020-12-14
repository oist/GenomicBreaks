# Service functions that change environments in an horrific un-sane way,
# but that I find useful at the command line.

pgeAndNxt <- function () {
  env <- environment()
  env$r <- NULL
  env$x <- NULL
  pge <- function(obj) {
    env$r <<- 1:10
    env$x <<- obj
    x[r]
  }
  nxt <- function () {
    env$r <<- r + 10
    x <- env$x[r]
    names(x) <- r
    x
  }
  list(pge,nxt)
}

tmpObj <- pgeAndNxt()
pge <- tmpObj[[1]]
nxt <- tmpObj[[2]]
rm(tmpObj)

Head <- function(...) {
  opt <- options("showTailLines" = Inf)
  print(head(...))
  options(opt) # Sometimes it does not restore the environment properly.
}
