% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strandNames.R
\name{strandNames}
\alias{strandNames}
\title{Sequence name with strand information}
\usage{
strandNames(gb, flip = FALSE, query = FALSE)
}
\arguments{
\item{gb}{A \code{GenomicBreaks} object of a \code{GRangesList} of \code{GenomicBreaks}
objects}
}
\value{
Returns a character vector, or a list of character vectors if the
input was a \code{GRangesList}.
}
\description{
Extract sequence names and pastes strand information to it.
}
\examples{
gb1       <- GenomicRanges::GRanges(c(A="Ref:100-200:+", B="Ref:400-500:+"))
gb1$query <- GenomicRanges::GRanges(c(A="Que:100-200",   B="Que:400-500"))
strandNames(gb1)
queryStrandNames(gb1)

}