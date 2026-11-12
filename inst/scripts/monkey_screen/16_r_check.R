suppressMessages(library(GenomicBreaks))
args <- commandArgs(trailingOnly = TRUE)
funs <- c("albert_index", "marcel_index", "rafiki_index", "wukong_index", "abu_index", "imo_index",
          "george_index", "jack_index", "baker_index", "hanuman_index", "mizaru_index", "able_index")
geom <- setdiff(funs, c("george_index", "baker_index"))
out <- file(args[1], "w")
for (f in readLines(args[2])) {
  t0 <- proc.time()[3]
  gb <- load_genomic_breaks(f)
  # filter of skeleton.Rmd
  if (length(seqlevels(gb)) > 250 || length(seqlevels(gb$query)) > 250) {
    gb <- gb[seqnames(gb) %in% names(head(sort(sum(split(width(gb), seqnames(gb))), decr=T), 250))]
    gb <- gb[seqnames(gb$query) %in% names(head(sort(sum(split(width(gb$query), seqnames(gb$query))), decr=T), 250))]
    seqlevels(gb) <- seqlevelsInUse(gb)
    seqlevels(gb$query) <- seqlevelsInUse(gb$query)
  }
  t1 <- proc.time()[3]
  v <- sapply(funs, function(fn) get(fn)(gb))
  t2 <- proc.time()[3]
  sw <- swap(gb, sort = TRUE)
  w <- sapply(geom, function(fn) get(fn)(sw))
  rv <- sapply(geom, function(fn) get(fn)(reverse(gb)))
  rq <- sapply(geom, function(fn) get(fn)(reverse(gb, query = TRUE)))
  cat(basename(f), length(gb), round(t1 - t0, 1), round(t2 - t1, 2), sprintf("%.15g", v), "|", sprintf("%.15g", w), "|",
      sprintf("%.15g", rv), "|", sprintf("%.15g", rq), "\n", file = out)
  flush(out)
}
close(out)
