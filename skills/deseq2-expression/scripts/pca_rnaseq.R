# Standard RNA-seq PCA from an RSEM gene expected_count matrix.
# Usage: Rscript pca_rnaseq.R <rsem_gene_expected_count.tsv> <samplesheet.tsv> <out_dir>
#   samplesheet.tsv : columns sample, group [, subset]; group order in the plots = order of first appearance in the samplesheet
# Method (DESeq2 convention, as in plotPCA): filter low genes
# (>= MIN_COUNT counts in >= MIN_SAMPLES samples) -> size factors + variance stabilising transformation (vst, blind = TRUE, so no design
# influences it) -> PCA (prcomp, centred, not scaled) on the NTOP most variable genes, and on all genes.
# Also: PC1-PC3 / PC2-PC3, scree plot, sample-to-sample distance heatmap, top loading genes per PC.
# env (optional): MIN_COUNT (10), MIN_SAMPLES (3), NTOP (500)
suppressPackageStartupMessages({library(DESeq2); library(ggplot2); library(ggrepel)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript pca_rnaseq.R <matrix.tsv> <samplesheet.tsv> <out_dir>")
mat_file <- args[1]; ss_file <- args[2]; out_dir <- args[3]
min_count <- as.integer(Sys.getenv("MIN_COUNT", "10")); min_samples <- as.integer(Sys.getenv("MIN_SAMPLES", "3")); ntop <- as.integer(Sys.getenv("NTOP", "500"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
m <- read.delim(mat_file, check.names = FALSE, stringsAsFactors = FALSE)
rownames(m) <- m$id; m$id <- NULL
ss <- read.delim(ss_file, stringsAsFactors = FALSE, colClasses = "character")
stopifnot(c("sample", "group") %in% names(ss), !anyDuplicated(ss$sample))
if (!"subset" %in% names(ss)) ss$subset <- "all"
miss <- setdiff(ss$sample, colnames(m)); if (length(miss)) stop("samples not in the matrix: ", paste(miss, collapse = ", "))
cts <- as.matrix(round(m[, ss$sample])); storage.mode(cts) <- "integer"
symbol <- ifelse(grepl("_", rownames(cts)), sub("^[^_]*_", "", rownames(cts)), rownames(cts)); names(symbol) <- rownames(cts)
coldata <- data.frame(sample = ss$sample, group = factor(ss$group, levels = unique(ss$group)), subset = factor(ss$subset, levels = unique(ss$subset)),
                      row.names = ss$sample)
keep <- rowSums(cts >= min_count) >= min_samples
cat(sum(keep), "of", nrow(cts), "genes kept (>=", min_count, "counts in >=", min_samples, "of", ncol(cts), "samples)\n")
dds <- DESeqDataSetFromMatrix(countData = cts[keep, ], colData = coldata, design = ~ 1)
dds <- estimateSizeFactors(dds)
write.csv(data.frame(sample = colnames(dds), size_factor = sizeFactors(dds), total_counts = colSums(counts(dds))), file.path(out_dir, "size_factors.csv"), row.names = FALSE)
vsd <- vst(dds, blind = TRUE)
X <- assay(vsd)
write.csv(X, file.path(out_dir, "vst_matrix.csv"))

grp <- levels(coldata$group)
pal <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
cols <- setNames(if (length(grp) <= length(pal)) pal[seq_along(grp)] else hcl.colors(length(grp), "Dark 3"), grp)
pca_run <- function(X, ntop, tag) {
  rv <- apply(X, 1, var)
  sel <- if (is.null(ntop)) rownames(X) else names(sort(rv, decreasing = TRUE))[seq_len(min(ntop, length(rv)))]
  p <- prcomp(t(X[sel, ]), center = TRUE, scale. = FALSE)
  pv <- 100 * p$sdev^2 / sum(p$sdev^2)
  sc <- data.frame(sample = rownames(p$x), p$x[, 1:min(6, ncol(p$x))], coldata[rownames(p$x), c("group", "subset")], row.names = NULL)
  write.csv(sc, file.path(out_dir, paste0("pca_scores_", tag, ".csv")), row.names = FALSE)
  write.csv(data.frame(PC = paste0("PC", seq_along(pv)), percent_variance = pv, cumulative = cumsum(pv)), file.path(out_dir, paste0("pca_variance_", tag, ".csv")), row.names = FALSE)
  ld <- p$rotation[, 1:min(4, ncol(p$rotation))]
  top <- do.call(rbind, lapply(colnames(ld), function(pc) { o <- order(abs(ld[, pc]), decreasing = TRUE)[1:min(25, nrow(ld))]
            data.frame(PC = pc, symbol = symbol[rownames(ld)[o]], loading = ld[o, pc], row.names = NULL) }))
  write.csv(top, file.path(out_dir, paste0("pca_top_loadings_", tag, ".csv")), row.names = FALSE)
  pc_plot <- function(a, b, fn) {
    g <- ggplot(sc, aes(.data[[a]], .data[[b]], color = group, shape = subset)) + geom_point(size = 3.4, alpha = 0.9) +
      geom_text_repel(aes(label = sample), size = 2.4, color = "grey25", max.overlaps = 100, min.segment.length = 0.2, show.legend = FALSE) +
      scale_color_manual(values = cols) + theme_minimal() +
      labs(title = paste0("PCA (", tag, ")"), x = sprintf("%s (%.1f%%)", a, pv[as.integer(sub("PC", "", a))]), y = sprintf("%s (%.1f%%)", b, pv[as.integer(sub("PC", "", b))]), color = "group", shape = "subset")
    ggsave(file.path(out_dir, paste0(fn, ".png")), g, width = 8, height = 6.2, dpi = 300); ggsave(file.path(out_dir, paste0(fn, ".pdf")), g, width = 8, height = 6.2)
  }
  pc_plot("PC1", "PC2", paste0("pca_PC1_PC2_", tag)); pc_plot("PC1", "PC3", paste0("pca_PC1_PC3_", tag)); pc_plot("PC2", "PC3", paste0("pca_PC2_PC3_", tag))
  npc <- min(10, length(pv))
  sv <- data.frame(PC = factor(paste0("PC", 1:npc), levels = paste0("PC", 1:npc)), pv = pv[1:npc])
  g <- ggplot(sv, aes(PC, pv)) + geom_col(fill = "#4c78a8", width = 0.7) + geom_text(aes(label = sprintf("%.1f", pv)), vjust = -0.4, size = 3) + theme_minimal() +
    labs(title = paste0("Variance explained (", tag, ")"), x = NULL, y = "% variance")
  ggsave(file.path(out_dir, paste0("pca_scree_", tag, ".png")), g, width = 6.5, height = 4.2, dpi = 300)
  cat(tag, ": genes used", length(sel), "; PC1-PC4 variance", paste(sprintf("%.1f%%", pv[1:4]), collapse = " "), "\n")
  invisible(list(sc = sc, pv = pv))
}
r_top <- pca_run(X, ntop, paste0("top", ntop))
r_all <- pca_run(X, NULL, "allgenes")

# sample-to-sample distance heatmap (Euclidean on vst, hierarchical order)
d <- dist(t(X)); hc <- hclust(d); ord <- hc$labels[hc$order]
dm <- as.matrix(d)
dl <- data.frame(a = factor(rep(rownames(dm), ncol(dm)), levels = ord), b = factor(rep(colnames(dm), each = nrow(dm)), levels = rev(ord)), dist = as.vector(dm))
g <- ggplot(dl, aes(a, b, fill = dist)) + geom_tile() + scale_fill_viridis_c(direction = -1, name = "Euclidean\ndistance (vst)") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), axis.text.y = element_text(size = 7)) + labs(title = "Sample-to-sample distances", x = NULL, y = NULL)
ggsave(file.path(out_dir, "sample_distance_heatmap.png"), g, width = 8.5, height = 7.5, dpi = 300); ggsave(file.path(out_dir, "sample_distance_heatmap.pdf"), g, width = 8.5, height = 7.5)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("DONE\n")
