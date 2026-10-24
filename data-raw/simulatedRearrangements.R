library("GenomicBreaks")

# A set of simulated genome pairs for testing and benchmarking the functions
# that flag structural variants.  Each simulation is stored with base R types
# only, so that the file does not depend on the serialisation format of
# Bioconductor objects.  Rebuild a GBreaks object with:
#
#   with(sims[[1]], permutationVectorToGBreaks(v, tsplit, qsplit))
#
# and recover the events from sims[[1]]$events.

set.seed(20261024)

nsim <- 200
sims <- lapply(seq_len(nsim), function(i) {
  params <- list(n           = 14,
                 nEvents     = sample(1:4, 1),
                 nTargetSeqs = sample(1:2, 1),
                 nQuerySeqs  = sample(1:2, 1))
  gb <- do.call(simulateRearrangements, params)
  list(v      = permutationVector(gb),
       tsplit = as.vector(table(seqnames(gb))),
       qsplit = as.vector(table(seqnames(gb$query))),
       events = S4Vectors::metadata(gb)$events,
       params = params)
})
names(sims) <- sprintf("sim%03d", seq_len(nsim))

saveRDS(sims, "inst/extdata/simulatedRearrangements.rds", compress = "xz")
