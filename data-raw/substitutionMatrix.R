exampleSubstitutionMatrix <- structure(c(6635437L, 414321L, 846130L, 512526L,
0L, 661258L, 424914L, 4887709L, 359307L, 909347L, 0L, 387445L, 911388L, 359509L,
4892634L, 423292L, 0L, 387195L, 513318L, 845040L, 413602L, 6626372L, 0L,
661137L, 0L, 0L, 0L, 0L, 0L, 0L, 851577L, 535002L, 534881L, 850529L, 0L, 0L),
dim = c(6L, 6L), dimnames = list(c("A", "C", "G", "T", "N", "-"), c("A", "C", "G", "T", "N", "-")))
usethis::use_data(exampleSubstitutionMatrix, overwrite = TRUE)

exampleProbabilityMatrix <-  readTrainFile(system.file("extdata/example.train", package = "GenomicBreaks")) $ probability_matrix
usethis::use_data(exampleProbabilityMatrix, overwrite = TRUE)
