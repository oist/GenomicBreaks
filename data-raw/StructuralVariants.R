library("GenomicBreaks")

# First, let's define colinear regions

exampleColinear                           <- GRanges(c("chrA:100-150", "chrA:251-300"))
strand(exampleColinear)                   <- c(              "+",            "+")
exampleColinear$query                     <- GRanges(c("chrB:100-150", "chrB:251-300"))
exampleColinear                           <- GBreaks(exampleColinear)
seqlengths(exampleColinear)               <- seqlengths(exampleColinear$query) <- 600
isSorted(exampleColinear)

# And a counter-example

exampleNotColinear                        <- GRanges(c("chrA:100-150", "chrA:251-300"))
strand(exampleNotColinear)                <- c(              "+",            "+")
exampleNotColinear$query                  <- GRanges(c("chrB:201-251",  "chrB:50-100"))
exampleNotColinear                        <- GBreaks(exampleNotColinear)
seqlengths(exampleNotColinear)            <- seqlengths(exampleNotColinear$query) <- 600
isSorted(exampleNotColinear)

# A colinear triplet for sanity checks

exampleColinear3                          <- GRanges(c("chrA:100-200", "chrA:201-300", "chrA:301-400"))
strand(exampleColinear3)                  <- c(              "+",            "+",            "+")
exampleColinear3$query                    <- GRanges(c("chrB:100-200", "chrB:201-300", "chrB:301-400"))
exampleColinear3                          <- GBreaks(exampleColinear3)
seqlengths(exampleColinear3)              <- seqlengths(exampleColinear3$query) <- 600
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
seqlengths(exampleInversion)              <- seqlengths(exampleInversion$query) <- 600
isSorted(exampleInversion)

# Inversion in a 5-bloc context with non-collapsed regions
# This is an intermediary step towards `exampleDoubleInversion1`

# ABC/ABC -> ABC/AbC

exampleInversion5uncollapsed              <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590"))
strand(exampleInversion5uncollapsed)      <- c(              "+",            "-",            "-",            "+",            "+")
exampleInversion5uncollapsed$query        <- GRanges(c("chrB:100-190", "chrB:300-390", "chrB:200-290", "chrB:400-490", "chrB:500-590"))
exampleInversion5uncollapsed              <- GBreaks(exampleInversion5uncollapsed)
seqlengths(exampleInversion5uncollapsed)  <- seqlengths(exampleInversion5uncollapsed$query) <- 600
isSorted(exampleInversion5uncollapsed)

# Double inversions
# ABC/ABC -> ABC/baC -> ABC/bcA => 3+, 1-, 2-

exampleDoubleInversion1                   <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590"))
strand(exampleDoubleInversion1)           <- c(              "+",            "+",            "-",            "-",            "+")
exampleDoubleInversion1$query             <- GRanges(c("chrB:100-190", "chrB:400-490", "chrB:200-290", "chrB:300-390", "chrB:500-590"))
exampleDoubleInversion1                   <- GBreaks(exampleDoubleInversion1)
seqlengths(exampleDoubleInversion1)       <- seqlengths(exampleDoubleInversion1$query) <- 600
isSorted(exampleDoubleInversion1)

# Double inversions
# ABC/ABC -> ABC/Acb -> ABC/Cab => 2-, 3-, 1+
exampleDoubleInversion2                   <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590"))
strand(exampleDoubleInversion2)           <- c(              "+",            "-",            "-",            "+",            "+")
exampleDoubleInversion2$query             <- GRanges(c("chrB:100-190", "chrB:300-390", "chrB:400-490", "chrB:200-290", "chrB:500-590"))
exampleDoubleInversion2                   <- GBreaks(exampleDoubleInversion2)
seqlengths(exampleDoubleInversion2)       <- seqlengths(exampleDoubleInversion2$query) <- 600
isSorted(exampleDoubleInversion2)

# Double inversions
# ABC/cba -> ABC/BCa -> ABC/BAc => 2+, 1+, 3-
exampleDoubleInversion3                   <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590"))
strand(exampleDoubleInversion3)           <- c(              "-",            "+",            "+",            "-",            "-")
exampleDoubleInversion3$query             <- GRanges(c("chrB:500-590", "chrB:300-390", "chrB:200-290", "chrB:400-490", "chrB:100-190"))
exampleDoubleInversion3                   <- GBreaks(exampleDoubleInversion3)
seqlengths(exampleDoubleInversion3)       <- seqlengths(exampleDoubleInversion3$query) <- 600
isSorted(exampleDoubleInversion3)

# Double inversions
# ABC/cba -> ABC/cAB -> ABC/aCB => 1-, 3+, 2+

exampleDoubleInversion4                   <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590"))
strand(exampleDoubleInversion4)           <- c(              "-",            "-",            "+",            "+",            "-")
exampleDoubleInversion4$query             <- GRanges(c("chrB:500-590", "chrB:200-290", "chrB:400-490", "chrB:300-390", "chrB:100-190"))
exampleDoubleInversion4                   <- GBreaks(exampleDoubleInversion4)
seqlengths(exampleDoubleInversion4)       <- seqlengths(exampleDoubleInversion4$query) <- 600
isSorted(exampleDoubleInversion4)

# Nested inversions
# ABCDE -> AdcbE -> AdCbE

exampleNestedInversions                   <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590"))
strand(exampleNestedInversions)           <- c(              "+",            "-",            "+",            "-",            "+")
exampleNestedInversions$query             <- GRanges(c("chrA:100-190", "chrA:400-490", "chrA:300-390", "chrA:200-290", "chrA:500-590"))
exampleNestedInversions                   <- GBreaks(exampleNestedInversions)
seqlengths(exampleNestedInversions)       <- seqlengths(exampleNestedInversions$query) <- 600
isSorted(exampleNestedInversions)

# Twin inversions
# ABCD -> AbCD -> AbcD

exampleTwinInversions                          <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490"))
strand(exampleTwinInversions)                  <- c(              "+",            "-",            "-",            "+"      )
exampleTwinInversions$query                    <- GRanges(c("chrB:100-190", "chrB:200-290", "chrB:300-390", "chrB:400-490"))
exampleTwinInversions                          <- GBreaks(exampleTwinInversions)
seqlengths(exampleTwinInversions)              <- seqlengths(exampleTwinInversions$query) <- 600
isSorted(exampleTwinInversions)

exampleNotTwinInversions                          <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490"))
strand(exampleNotTwinInversions)                  <- c(              "+",            "-",            "-",            "-"      )
exampleNotTwinInversions$query                    <- GRanges(c("chrB:100-190", "chrB:200-290", "chrB:300-390", "chrB:400-490"))
exampleNotTwinInversions                          <- GBreaks(exampleNotTwinInversions)
seqlengths(exampleNotTwinInversions)              <- seqlengths(exampleNotTwinInversions$query) <- 600
isSorted(exampleNotTwinInversions)

# Translocation

exampleTranslocation                      <- GRanges(c("chrA:100-200", "chrA:201-300", "chrA:301-400"))
strand(exampleTranslocation)              <- c(              "+",            "+",            "+")
exampleTranslocation$query                <- GRanges(c("chrB:100-200", "chrC:201-300", "chrB:301-400"))
exampleTranslocation                      <- GBreaks(exampleTranslocation)
seqlengths(exampleTranslocation)          <- 600
seqlengths(exampleTranslocation$query)    <- c(600,600)
isSorted(exampleTranslocation)

# Translocation with minus strand

exampleTranslocation2                     <- GRanges(c("chrA:100-200", "chrA:201-300", "chrA:301-400"))
strand(exampleTranslocation2)             <- c(              "+",            "-",            "+")
exampleTranslocation2$query               <- GRanges(c("chrB:100-200", "chrC:201-300", "chrB:301-400"))
exampleTranslocation2                     <- GBreaks(exampleTranslocation2)
seqlengths(exampleTranslocation2)         <- 600
seqlengths(exampleTranslocation2$query)   <- c(600,600)
isSorted(exampleTranslocation2)

# Translocation downstream

exampleTranslocation3                     <- GRanges(c("chrA:100-200", "chrA:201-300", "chrA:301-400"))
strand(exampleTranslocation3)             <- c(              "+",            "+",            "+")
exampleTranslocation3$query               <- GRanges(c("chrB:100-200", "chrB:301-400", "chrB:201-300"))
exampleTranslocation3                     <- GBreaks(exampleTranslocation3)
seqlengths(exampleTranslocation3)         <- seqlengths(exampleTranslocation3$query) <- 600
isSorted(exampleTranslocation3)

# Insertion on the target

exampleInsertion                          <- GRanges(c("chrA:100-200", "chrC:401-500", "chrA:201-300"))
strand(exampleInsertion)                  <- c(              "+",            "+",            "+")
exampleInsertion$query                    <- GRanges(c("chrB:100-200", "chrB:201-300", "chrB:301-400"))
exampleInsertion                          <- GBreaks(exampleInsertion)
seqlengths(exampleInsertion)              <- c(600,600)
seqlengths(exampleInsertion$query)        <- 600
exampleInsertion                          <- sort(exampleInsertion, ignore.strand = TRUE)
isSorted(exampleInsertion)

# Deletion on the target

exampleDeletion                           <- GRanges(c("chrA:100-200", "chrA:201-300", "chrA:301-400"))
strand(exampleDeletion)                   <- c(              "+",            "+",            "+")
exampleDeletion$query                     <- GRanges(c("chrB:100-200", "chrC:401-500", "chrB:201-300"))
exampleDeletion                           <- GBreaks(exampleDeletion)
seqlengths(exampleDeletion)               <- 600
seqlengths(exampleDeletion$query)         <- c(600,600)
isSorted(exampleDeletion)


usethis::use_data(
  exampleColinear,
  exampleNotColinear,
  exampleColinear3,
  exampleColinear5,
  exampleInversion,
  exampleInversion5uncollapsed,
  exampleDoubleInversion1,
  exampleDoubleInversion2,
  exampleDoubleInversion3,
  exampleDoubleInversion4,
  exampleNestedInversions,
  exampleTwinInversions,
  exampleNotTwinInversions,
  exampleTranslocation,
  exampleTranslocation2,
  exampleInsertion,
  exampleDeletion,
  overwrite = TRUE)
