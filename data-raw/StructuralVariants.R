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

# Inversions 2

exampleInversion2                          <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590", "chrA:600-690"))
strand(exampleInversion2)                  <- c(              "+",            "-",            "-",            "-",            "+",            "+"      )
exampleInversion2$query                    <- GRanges(c("chrA:100-190", "chrA:300-390", "chrA:200-290", "chrA:400-490", "chrA:600-690", "chrA:500-590"))
exampleInversion2                          <- GBreaks(exampleInversion2)
seqlengths(exampleInversion2)              <- seqlengths(exampleInversion2$query) <- 700
isSorted(exampleInversion2)

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

# Example of inversions from literature.

# Example used in the paper from Garg et al. (2019).
# {-2,5,4,-1,3,6,9,-7,-8}
# Reversal distance = 5 reversals. 
# Details for Reversal distance (d) computation: 10 breakpoints (b); 5 cycles(c); 0 hurdles(h): d = b-c+h (+1 if fortress).
exampleInversionGarg2019              <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590", "chrA:600-690", "chrA:700-790", "chrA:800-890", "chrA:900-990"))
strand(exampleInversionGarg2019)      <- c(              "-",            "+",            "+",            "-",            "+",            "+",            "+",            "-",            "-")
exampleInversionGarg2019$query        <- GRanges(c("chrB:200-290", "chrB:500-590", "chrB:400-490", "chrB:100-190", "chrB:300-390", "chrB:600-690", "chrB:900-990", "chrB:700-790", "chrB:800-890"))
exampleInversionGarg2019              <- GBreaks(exampleInversionGarg2019)
seqlengths(exampleInversionGarg2019)  <- seqlengths(exampleInversionGarg2019$query) <- 1000
isSorted(exampleInversionGarg2019)


# Example used in the paper from Bader et al. (2001).
# (+3, +9, −7, +5, −10, +8, +4, −6, +11, +2, +1)
# - Reversal distance (d) = expected: 7 reversals.
exampleInversionBader2001             <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590",   "chrA:600-690", "chrA:700-790", "chrA:800-890", "chrA:900-990",   "chrA:1000-1090", "chrA:1100-1190"))
strand(exampleInversionBader2001)     <- c(              "+",            "+",            "-",            "+",            "-",            "+",            "+",            "-",            "+",             "+",             "+")
exampleInversionBader2001$query       <- GRanges(c("chrB:300-390", "chrB:900-990", "chrB:700-790", "chrB:500-590", "chrB:1000-1090", "chrB:800-890", "chrB:400-490", "chrB:600-690", "chrB:1100-1190", "chrB:200-290",   "chrB:100-190"))
exampleInversionBader2001             <- GBreaks(exampleInversionBader2001)
seqlengths(exampleInversionBader2001) <- seqlengths(exampleInversionBader2001$query) <- 1200
isSorted(exampleInversionBader2001)


# Example used in the paper from Hannehalli and Pevzner (1999) (Figure 4(a)).
# {+5, +7, +6, +8, +1, +3, +2, +4}
# - Reversal distance (d) = expected: 8 reversals
# Details for Reversal distance (d) computation: 9 breakpoints (b); 3 cycles(c); 2 hurdles(h): d = b-c+h (+1 if fortress).
exampleInversionHP1999fig4a             <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590", "chrA:600-690", "chrA:700-790", "chrA:800-890"))
strand(exampleInversionHP1999fig4a)     <- c(              "+",            "+",            "+",            "+",            "+",            "+",            "+",            "+")
exampleInversionHP1999fig4a$query       <- GRanges(c("chrB:500-590", "chrB:700-790", "chrB:600-690", "chrB:800-890", "chrB:100-190", "chrB:300-390", "chrB:200-290", "chrB:400-490"))
exampleInversionHP1999fig4a             <- GBreaks(exampleInversionHP1999fig4a)
seqlengths(exampleInversionHP1999fig4a) <- seqlengths(exampleInversionHP1999fig4a$query) <- 900
isSorted(exampleInversionHP1999fig4a)


# Example used in the paper from Hannehalli and Pevzner (1999) (Figure 4(b)).
# {+2, +4, +3, +5, +7, +6, +8, +1}
# - Reversal distance (d) = expected: 9 reversals
# Details for Reversal distance (d) computation: 9 breakpoints (b); 3 cycles(c); 3 hurdles(h): d = b-c+h (+1 if fortress).
exampleInversionHP1999fig4b             <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590", "chrA:600-690", "chrA:700-790", "chrA:800-890"))
strand(exampleInversionHP1999fig4b)     <- c(              "+",            "+",            "+",            "+",            "+",            "+",            "+",            "+")
exampleInversionHP1999fig4b$query       <- GRanges(c("chrB:200-290", "chrB:400-490", "chrB:300-390", "chrB:500-590", "chrB:700-790", "chrB:600-690", "chrB:800-890", "chrB:100-190"))
exampleInversionHP1999fig4b             <- GBreaks(exampleInversionHP1999fig4b)
seqlengths(exampleInversionHP1999fig4b) <- seqlengths(exampleInversionHP1999fig4b$query) <- 900
isSorted(exampleInversionHP1999fig4b)

# Example used in the book ``Mathematics of Evolution and Phylogeny`` (2005) (Section 10.4.2).
# {0, 2, 1, 3, 5, 7, 6, 8, 9, 4, 10}
# - Reversal distance (d) = expected: 8 reversals
# Details for Reversal distance (d) computation: 9 breakpoints (b); 4 cycles(c); 2 hurdles(h): d = b-c+h (+1 if fortress).
exampleInversionBergeron2005a             <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590", "chrA:600-690", "chrA:700-790", "chrA:800-890", "chrA:900-990",   "chrA:1000-1090", "chrA:1100-1190"))
strand(exampleInversionBergeron2005a)     <- c(              "+",            "+",            "+",            "+",            "+",            "+",            "+",            "+",            "+",             "+",             "+")
exampleInversionBergeron2005a$query       <- GRanges(c("chrB:100-190", "chrB:300-390", "chrB:200-290", "chrB:400-490", "chrB:600-690", "chrB:800-890", "chrB:700-790", "chrB:900-990", "chrB:1000-1090", "chrB:500-590",   "chrB:1100-1190"))
exampleInversionBergeron2005a             <- GBreaks(exampleInversionBergeron2005a)
seqlengths(exampleInversionBergeron2005a) <- seqlengths(exampleInversionBergeron2005a$query) <- 1200
isSorted(exampleInversionBergeron2005a)


# Example used in the book ``Mathematics of Evolution and Phylogeny`` (2005) (Figure 10.6).
# P_2 = {0, -3, 1, 2, 4, 6, 5, 7, -15, -13, -14, -12, -10, -11, -9, 8, 16}
# - Reversal distance (d) = expected: 13 reversals
# Details for Reversal distance (d) computation: 15 breakpoints (b); 6 cycles(c); 3 hurdles(h): d = b-c+h (+1 if fortress).
exampleInversionBergeron2005b             <- GRanges(c("chrA:100-190", "chrA:200-290", "chrA:300-390", "chrA:400-490", "chrA:500-590", "chrA:600-690", "chrA:700-790", "chrA:800-890", "chrA:900-990",   "chrA:1000-1090", "chrA:1100-1190",  "chrA:1200-1290", "chrA:1300-1390", "chrA:1400-1490", "chrA:1500-1590", "chrA:1600-1690", "chrA:1700-1790"))
strand(exampleInversionBergeron2005b)     <- c(              "+",            "-",           "+",            "+",            "+",            "+",            "+",            "+",              "-",              "-",             "-",               "-",             "-",               "-",             "-",              "+",              "+")
exampleInversionBergeron2005b$query       <- GRanges(c("chrB:100-190", "chrB:400-490", "chrB:200-290", "chrB:300-390", "chrB:500-590", "chrB:700-790", "chrB:600-690", "chrB:800-890", "chrB:1600-1690", "chrB:1400-1490",  "chrB:1500-1590", "chrB:1300-1390", "chrB:1100-1190", "chrB:1200-1290", "chrB:1000-1090", "chrB:900-990",   "chrB:1700-1790"))
exampleInversionBergeron2005b             <- GBreaks(exampleInversionBergeron2005b)
seqlengths(exampleInversionBergeron2005b) <- seqlengths(exampleInversionBergeron2005b$query) <- 1800
isSorted(exampleInversionBergeron2005b)

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
  exampleInversion2,
  exampleInversion5uncollapsed,
  exampleDoubleInversion1,
  exampleDoubleInversion2,
  exampleDoubleInversion3,
  exampleDoubleInversion4,
  exampleNestedInversions,
  exampleTwinInversions,
  exampleNotTwinInversions,
  exampleInversionGarg2019,
  exampleInversionBader2001,
  exampleInversionHP1999fig4a,
  exampleInversionHP1999fig4b,
  exampleInversionBergeron2005a,
  exampleInversionBergeron2005b,
  exampleTranslocation,
  exampleTranslocation2,
  exampleInsertion,
  exampleDeletion,
  overwrite = TRUE)
