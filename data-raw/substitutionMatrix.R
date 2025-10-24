exampleSubstitutionMatrix <-  readTrainFile(system.file("extdata/example.train", package = "GenomicBreaks")) $ probability_matrix
usethis::use_data(exampleSubstitutionMatrix, overwrite = TRUE)
