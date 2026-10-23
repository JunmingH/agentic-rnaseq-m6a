# Guitar transcript-wide distribution (5'UTR / CDS / 3'UTR metagene) of DART-seq sites, RAC-filter and no-RAC-filter versions.
# Usage: Rscript 02_guitar.R <gtf> <bed_dir> <out_dir> [organism]
#   bed_dir: output of 01_make_beds.py (<group>_{norac,rac}.bed, groups.txt = group order in the plot); organism default "Homo sapiens"
# GuitarPlot settings as in our earlier DART-seq distribution plots ( headOrtail = TRUE, enableCI = FALSE, mapFilterTranscript = TRUE, pltTxType = "mrna"); the TxDb is
#   built once from the SAME gtf as the STAR index and passed to both calls. One call per version with all groups -> overlay plot, plotted densities, plot object.
suppressPackageStartupMessages({library(GenomicFeatures); library(Guitar); library(ggplot2)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript 02_guitar.R <gtf> <bed_dir> <out_dir> [organism]")
gtf <- args[1]; bed_dir <- args[2]; out_dir <- args[3]; organism <- if (length(args) >= 4) args[4] else "Homo sapiens"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE); setwd(out_dir)
GROUPS <- readLines(file.path(bed_dir, "groups.txt")); GROUPS <- GROUPS[GROUPS != ""]
cat("groups:", paste(GROUPS, collapse = ", "), "\nbuilding TxDb from", gtf, "\n")
txdb <- makeTxDbFromGFF(gtf, format = "gtf", organism = organism)
# (txGuitarTxdb expects a file that Guitar reads with read.table, not the object from makeGuitarTxdb; so the TxDb object is passed and
#  GuitarPlot builds the Guitar coordinates itself)
for (ver in c("norac", "rac")) {
  beds <- file.path(bed_dir, paste0(GROUPS, "_", ver, ".bed")); stopifnot(all(file.exists(beds)))
  cat("\n=== GuitarPlot", ver, "\n")
  # (no miscOutFilePrefix: with it GuitarPlot only writes a pdf and returns NULL)
  p <- GuitarPlot(txTxdb = txdb, stBedFiles = beds, stGroupName = GROUPS,
                  headOrtail = TRUE, enableCI = FALSE, mapFilterTranscript = TRUE, pltTxType = "mrna")
  stopifnot(inherits(p, "ggplot"))
  ggsave(file.path(out_dir, paste0("guitar_", ver, "_guitar_default.pdf")), p, width = 8, height = 6)   # Guitar's own look and size
  d <- p$data
  write.table(d, file.path(out_dir, paste0("guitar_", ver, "_plot_data.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  ttl <- paste0("DART-seq sites on mRNA (", if (ver == "rac") "RAC filter" else "no RAC filter", ")")
  p1 <- p + theme_bw() + labs(title = ttl)
  ggsave(file.path(out_dir, paste0("guitar_", ver, "_overlay.png")), p1, width = 9, height = 5, dpi = 300)
  ggsave(file.path(out_dir, paste0("guitar_", ver, "_overlay.pdf")), p1, width = 9, height = 5)
  # Guitar's layers refer to its own data vectors, so neither facets nor %+% work on this object; plot subsets from *_plot_data.tsv instead.
  saveRDS(p, file.path(out_dir, paste0("guitar_", ver, "_plot.rds")))
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("DONE\n")
