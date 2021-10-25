library("GenomicBreaks")

# First, let's define colinear regions

exampleColinear                           <- GRanges(c("chrA:100-150", "chrA:251-300"))
strand(exampleColinear)                   <- c(              "+",            "+")
exampleColinear$query                     <- GRanges(c("chrB:100-150", "chrB:251-300"))
exampleColinear                           <- GBreaks(exampleColinear)
seqlengths(exampleColinear)               <- seqlengths(exampleColinear$query) <- 300
isSorted(exampleColinear)

# And a counter-example

exampleNotColinear                        <- GRanges(c("chrA:100-150", "chrA:251-300"))
strand(exampleNotColinear)                <- c(              "+",            "+")
exampleNotColinear$query                  <- GRanges(c("chrB:201-251",  "chrB:50-100"))
exampleNotColinear                        <- GBreaks(exampleNotColinear)
seqlengths(exampleNotColinear)            <- seqlengths(exampleNotColinear$query) <- 300
isSorted(exampleNotColinear)

# A colinear triplet for sanity checks

exampleColinear3                          <- GRanges(c("chrA:100-200", "chrA:201-300", "chrA:301-400"))
strand(exampleColinear3)                  <- c(              "+",            "+",            "+")
exampleColinear3$query                    <- GRanges(c("chrB:100-200", "chrB:201-300", "chrB:301-400"))
exampleColinear3                          <- GBreaks(exampleColinear3)
seqlengths(exampleColinear3)              <- seqlengths(exampleColinear3$query) <- 500
isSorted(exampleColinear3)

# A colinear 5-uplet for later example of double inversion

exampleColinear5                          <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590"))
strand(exampleColinear5)                  <- c(              "+",            "+",            "+",            "+",            "+")
exampleColinear5$query                    <- GRanges(c("chrB:100-190", "chrB:200-290", "chrB:300-390", "chrB:400-490", "chrB:500-590"))
exampleColinear5                          <- GBreaks(exampleColinear5)
seqlengths(exampleColinear5)              <- seqlengths(exampleColinear5$query) <- 600
isSorted(exampleColinear5)

# Inversions

exampleInversion                          <- GRanges(c("chrA:100-190", "chrA:210-291", "chrA:301-400"))
strand(exampleInversion)                  <- c(              "+",            "-",            "+")
exampleInversion$query                    <- GRanges(c("chrB:100-190", "chrB:210-291", "chrB:301-400"))
exampleInversion                          <- GBreaks(exampleInversion)
seqlengths(exampleInversion)              <- seqlengths(exampleInversion$query) <- 500
isSorted(exampleInversion)

# Inversion in a 5-bloc context with non-collapsed regions
# This is an intermediary step towards `exampleDoubleInversion`

exampleInversion5uncollapsed              <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590"))
strand(exampleInversion5uncollapsed)      <- c(              "+",            "-",            "-",            "+",            "+")
exampleInversion5uncollapsed$query        <- GRanges(c("chrB:100-190", "chrB:300-390", "chrB:200-290", "chrB:400-490", "chrB:500-590"))
exampleInversion5uncollapsed              <- GBreaks(exampleInversion5uncollapsed)
seqlengths(exampleInversion5uncollapsed)  <- seqlengths(exampleInversion5uncollapsed$query) <- 600
isSorted(exampleInversion5uncollapsed)

# Double inversions

exampleDoubleInversion                    <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590"))
strand(exampleDoubleInversion)            <- c(              "+",            "+",            "-",            "-",            "+")
exampleDoubleInversion$query              <- GRanges(c("chrB:100-190", "chrB:400-490", "chrB:200-290", "chrB:300-390", "chrB:500-590"))
exampleDoubleInversion                    <- GBreaks(exampleDoubleInversion)
seqlengths(exampleDoubleInversion)        <- seqlengths(exampleDoubleInversion$query) <- 600
isSorted(exampleDoubleInversion)

# Translocation

exampleTranslocation                      <- GRanges(c("chrA:100-200", "chrA:201-300", "chrA:301-400"))
strand(exampleTranslocation)              <- c(              "+",            "+",            "+")
exampleTranslocation$query                <- GRanges(c("chrB:100-200", "chrC:201-300", "chrB:301-400"))
exampleTranslocation                      <- GBreaks(exampleTranslocation)
seqlengths(exampleTranslocation)          <- 500
seqlengths(exampleTranslocation$query)    <- c(500,500)
isSorted(exampleTranslocation)

# Translocation with minus strand

exampleTranslocation2                     <- GRanges(c("chrA:100-200", "chrA:201-300", "chrA:301-400"))
strand(exampleTranslocation2)             <- c(              "+",            "-",            "+")
exampleTranslocation2$query               <- GRanges(c("chrB:100-200", "chrC:201-300", "chrB:301-400"))
exampleTranslocation2                     <- GBreaks(exampleTranslocation2)
seqlengths(exampleTranslocation2)         <- 500
seqlengths(exampleTranslocation2$query)   <- c(500,500)
isSorted(exampleTranslocation2)

# Insertion on the target

exampleInsertion                          <- GRanges(c("chrA:100-200", "chrC:401-500", "chrA:201-300"))
strand(exampleInsertion)                  <- c(              "+",            "+",            "+")
exampleInsertion$query                    <- GRanges(c("chrB:100-200", "chrB:201-300", "chrB:301-400"))
exampleInsertion                          <- GBreaks(exampleInsertion)
seqlengths(exampleInsertion)              <- c(500,500)
seqlengths(exampleInsertion$query)        <- 500
exampleInsertion                          <- sort(exampleInsertion, ignore.strand = TRUE)
isSorted(exampleInsertion)

# Deletion on the target

exampleDeletion                           <- GRanges(c("chrA:100-200", "chrA:201-300", "chrA:301-400"))
strand(exampleDeletion)                   <- c(              "+",            "+",            "+")
exampleDeletion$query                     <- GRanges(c("chrB:100-200", "chrC:401-500", "chrB:201-300"))
exampleDeletion                           <- GBreaks(exampleDeletion)
seqlengths(exampleDeletion)               <- 500
seqlengths(exampleDeletion$query)         <- c(500,500)
isSorted(exampleDeletion)


usethis::use_data(
  exampleColinear,
  exampleNotColinear,
  exampleColinear3,
  exampleColinear5,
  exampleInversion,
  exampleInversion5uncollapsed,
  exampleDoubleInversion,
  exampleTranslocation,
  exampleTranslocation2,
  exampleInsertion,
  exampleDeletion,
  overwrite = TRUE)
