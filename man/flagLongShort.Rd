% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flagLongShort.R
\name{flagLongShort}
\alias{flagLongShort}
\title{Flag long and short arms}
\usage{
flagLongShort(gr, annot, select = c("first", "last", "arbitrary"))
}
\arguments{
\item{gr}{A \code{\link{GBreaks}} or a simple \code{\link[GenomicRanges:GRanges-class]{GenomicRanges::GRanges}} object}

\item{annot}{A \code{GRanges} file containing the coordinate of arms and their
nature (such as \code{long}, \code{short}, \code{XSR} or \code{YSR}) in a \code{Type}
\emph{metadata column}.}

\item{select}{One method among \code{first}, \code{last} and \code{arbitrary}, to decide
on how to handle the regions that match both arms.}
}
\value{
Returns a modified version of the object with an extra \emph{metadata
column}, \code{Arm} in which the information from the annotation file's \code{Type}
column was transferred.  See the examples below and in the manual of
\code{\link[IRanges:findOverlaps-methods]{IRanges::findOverlaps}} for details. on how regions that match both arms
are handled.
}
\description{
\emph{Oikopleura} genomes are special as the long and short arms of their
chromosomes have different properties such as \code{GC} or repeat content.  It
can be useful to know if genomic region belongs to long or a short arm.
}
\examples{
annot <- GRanges(c("chrA:1-140", "chrA:150-500"), Type = c("short", "long"))

flagLongShort(exampleColinear, annot)
flagLongShort(exampleColinear, annot, select = "last")
flagLongShort(exampleColinear, annot, select = "arbitrary")

}
\seealso{
Other Flagging functions: 
\code{\link{flagAll}()},
\code{\link{flagColinearAlignments}()},
\code{\link{flagDoubleInversions}()},
\code{\link{flagInversions}()},
\code{\link{flagPairs}()},
\code{\link{flagTranslocations}()}
}
\author{
Charles Plessy
}
\concept{Flagging functions}
