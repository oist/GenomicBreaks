% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Translocations.R
\name{showTranslocations}
\alias{showTranslocations}
\title{Show translocations and their flanking blocks.}
\usage{
showTranslocations(gb, rename = TRUE)
}
\arguments{
\item{gb}{A \code{\link{GBreaks}} object processed with \code{\link{flagTranslocations}}.}

\item{rename}{Replace range names by their numeric order before subsetting.}
}
\value{
Returns the \code{GBreaks} object in which all ranges that are not
part of a translocation triplet have been discarded.  If the object was missing
the \code{tra} metadata column, return the object after discarding all of its
ranges.
}
\description{
Show translocations and their flanking blocks.
}
\examples{
showTranslocations(flagTranslocations(exampleTranslocation))

}
\seealso{
Other Translocation functions: 
\code{\link{filterTranslocations}()},
\code{\link{flagTranslocations}()}
}
\concept{Translocation functions}
