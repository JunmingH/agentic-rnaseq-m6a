# Differential expression (DESeq2) from an RSEM gene expected_count matrix (rsem-quantification skill, make_matrix.py).
# Usage: Rscript diff_expr.R <rsem_gene_expected_count.tsv> <samplesheet.tsv> <comparisons.tsv> <out_dir>
#   samplesheet.tsv : columns sample, group [, subset]   (sample = column name in the matrix; subset = e.g. cell line, optional)
#   comparisons.tsv : columns name, treated, control [, subset, label]
#                     treated/control = group names; subset = only the samples of that subset are modelled (empty or "all" = all samples)
# Per comparison: DESeqDataSetFromMatrix(design = ~ group) on the
# samples of its subset, control group = reference level, DESeq(betaPrior = FALSE), results(name = group_<treated>_vs_<control>, alpha),
# lfcShrink(apeglm) as an extra column, volcano (ggplot2/ggrepel). One fit is shared by comparisons with the same subset and control.
# Gene filter per fit: >= MIN_COUNT counts in at least as many samples as the smallest group of the fit.
# env (optional): ALPHA (0.05), LFC (1, volcano / up-down counts only), MIN_COUNT (10)
suppressPackageStartupMessages({library(DESeq2); library(ggplot2); library(ggrepel); library(openxlsx)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Usage: Rscript diff_expr.R <matrix.tsv> <samplesheet.tsv> <comparisons.tsv> <out_dir>")
mat_file <- args[1]; ss_file <- args[2]; cmp_file <- args[3]; out_dir <- args[4]
for (d in c("tables", "volcano", "qc")) dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
alpha <- as.numeric(Sys.getenv("ALPHA", "0.05")); lfc_cut <- as.numeric(Sys.getenv("LFC", "1")); min_count <- as.integer(Sys.getenv("MIN_COUNT", "10"))
cat("ALPHA", alpha, " LFC", lfc_cut, " MIN_COUNT", min_count, "\n")

m <- read.delim(mat_file, check.names = FALSE, stringsAsFactors = FALSE)
rownames(m) <- m$id; m$id <- NULL
ss <- read.delim(ss_file, stringsAsFactors = FALSE, colClasses = "character")
stopifnot(c("sample", "group") %in% names(ss), !anyDuplicated(ss$sample))
if (!"subset" %in% names(ss)) ss$subset <- "all"
miss <- setdiff(ss$sample, colnames(m)); if (length(miss)) stop("samples not in the matrix: ", paste(miss, collapse = ", "))
cts <- as.matrix(round(m[, ss$sample])); storage.mode(cts) <- "integer"      # RSEM expected_count is not an integer
symbol <- ifelse(grepl("_", rownames(cts)), sub("^[^_]*_", "", rownames(cts)), rownames(cts))   # --append-names: ENSG..._SYMBOL
names(symbol) <- rownames(cts)
coldata <- data.frame(sample = ss$sample, group = ss$group, subset = ss$subset, row.names = ss$sample, stringsAsFactors = FALSE)
cat("matrix:", nrow(cts), "genes x", ncol(cts), "samples\n"); print(table(coldata$group))

cmps <- read.delim(cmp_file, stringsAsFactors = FALSE, colClasses = "character")
stopifnot(c("name", "treated", "control") %in% names(cmps), !anyDuplicated(cmps$name))
if (!"subset" %in% names(cmps)) cmps$subset <- "all"
cmps$subset[is.na(cmps$subset) | cmps$subset == ""] <- "all"
if (!"label" %in% names(cmps)) cmps$label <- ""
cmps$label <- ifelse(is.na(cmps$label) | cmps$label == "", paste(cmps$treated, "vs", cmps$control), cmps$label)
for (g in c(cmps$treated, cmps$control)) if (!g %in% coldata$group) stop("comparison group not in the samplesheet: ", g)

# ---- QC: PCA of all samples (vst, blind)
keep_all <- rowSums(cts >= min_count) >= 3
dds_all <- DESeqDataSetFromMatrix(countData = cts[keep_all, ], colData = coldata, design = ~ 1)
vsd <- vst(dds_all, blind = TRUE)
pca <- prcomp(t(assay(vsd)))
pv <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pcd <- data.frame(pca$x[, 1:2], coldata[rownames(pca$x), c("group", "subset")], sample = rownames(pca$x))
write.csv(pcd, file.path(out_dir, "qc", "pca_coordinates.csv"), row.names = FALSE)
gp <- ggplot(pcd, aes(PC1, PC2, color = group, shape = subset)) + geom_point(size = 3.2) +
  geom_text_repel(aes(label = sample), size = 2.5, show.legend = FALSE, max.overlaps = 100) +
  theme_minimal() + labs(title = paste0("PCA of all ", ncol(cts), " samples (vst)"), x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"))
ggsave(file.path(out_dir, "qc", "pca_all_samples.png"), gp, width = 7.5, height = 6, dpi = 300)
ggsave(file.path(out_dir, "qc", "pca_all_samples.pdf"), gp, width = 7.5, height = 6)

# ---- comparisons
fits <- list(); summ <- data.frame(); sheets <- list()
for (i in seq_len(nrow(cmps))) {
  cp <- as.list(cmps[i, ])
  key <- paste(cp$subset, cp$control, sep = "|")
  if (is.null(fits[[key]])) {
    cd <- if (cp$subset == "all") coldata else coldata[coldata$subset == cp$subset, ]
    stopifnot(cp$control %in% cd$group)
    cd$group <- relevel(factor(cd$group), ref = cp$control)
    cc <- cts[, rownames(cd)]
    gsz <- min(table(cd$group))
    keep <- rowSums(cc >= min_count) >= gsz                # >= MIN_COUNT counts in at least as many samples as the smallest group
    cat("\n== fit", key, ":", ncol(cc), "samples, groups:", paste(levels(cd$group), collapse = "/"), "; genes kept", sum(keep), "of", nrow(cc), "\n")
    dds <- DESeqDataSetFromMatrix(countData = cc[keep, ], colData = cd, design = ~ group)
    dds <- DESeq(dds, betaPrior = FALSE)
    fits[[key]] <- dds
  }
  dds <- fits[[key]]
  coef <- paste0("group_", make.names(cp$treated), "_vs_", make.names(cp$control))
  if (!coef %in% resultsNames(dds)) stop("coefficient ", coef, " not in ", paste(resultsNames(dds), collapse = ", "))
  res <- results(dds, name = coef, alpha = alpha)
  shr <- lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)
  nc <- counts(dds, normalized = TRUE)
  g_num <- rownames(colData(dds))[colData(dds)$group == cp$treated]; g_den <- rownames(colData(dds))[colData(dds)$group == cp$control]
  tb <- data.frame(gene_id = sub("_.*$", "", rownames(res)), symbol = symbol[rownames(res)], full_id = rownames(res),
                   baseMean = res$baseMean, log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
                   log2FC_shrunk_apeglm = shr$log2FoldChange, pvalue = res$pvalue, padj = res$padj, stringsAsFactors = FALSE)
  tb[[paste0("mean_norm_", cp$treated)]] <- rowMeans(nc[, g_num, drop = FALSE])
  tb[[paste0("mean_norm_", cp$control)]] <- rowMeans(nc[, g_den, drop = FALSE])
  tb <- tb[order(tb$padj, na.last = TRUE), ]
  write.csv(tb, file.path(out_dir, "tables", paste0(cp$name, ".csv")), row.names = FALSE)
  sheets[[substr(cp$name, 1, 31)]] <- tb
  # ---- volcano
  tb$nlp <- -log10(tb$padj)
  ns_lab <- paste0("padj<", alpha, ", |log2FC|<", lfc_cut)
  tb$cls <- ifelse(is.na(tb$padj), "NA", ifelse(tb$padj < alpha & tb$log2FoldChange >= lfc_cut, "Up", ifelse(tb$padj < alpha & tb$log2FoldChange <= -lfc_cut, "Down",
                   ifelse(tb$padj < alpha, ns_lab, "n.s."))))
  pl <- tb[!is.na(tb$padj), ]
  pl$cls <- factor(pl$cls, levels = c("Up", "Down", ns_lab, "n.s."))
  lab <- pl[pl$padj < alpha & abs(pl$log2FoldChange) >= lfc_cut, ]
  lab <- head(lab[order(lab$padj), ], 20)
  tn <- table(pl$cls)
  vp <- ggplot(pl, aes(log2FoldChange, nlp, color = cls)) + geom_point(alpha = 0.7, size = 1.1) +
    scale_color_manual(values = setNames(c("#c0392b", "#2874a6", "#7f8c8d", "#d5d8dc"), c("Up", "Down", ns_lab, "n.s.")), drop = FALSE) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey40") +
    geom_text_repel(data = lab, aes(label = symbol), size = 2.6, color = "black", min.segment.length = 0, max.overlaps = 100, show.legend = FALSE) +
    theme_minimal() +
    labs(title = paste0("Volcano: ", cp$label), subtitle = paste0("DESeq2, n = ", length(g_num), " vs ", length(g_den), "; ", nrow(pl), " genes tested\nup ", tn[["Up"]], ", down ", tn[["Down"]], " (padj < ", alpha, " and |log2FC| >= ", lfc_cut, ")"),
         x = paste0("log2 fold change (", cp$treated, " / ", cp$control, ")"), y = "-log10(adjusted p)", color = NULL)
  ggsave(file.path(out_dir, "volcano", paste0("volcano_", cp$name, ".png")), vp, width = 7.5, height = 6.2, dpi = 300)
  ggsave(file.path(out_dir, "volcano", paste0("volcano_", cp$name, ".pdf")), vp, width = 7.5, height = 6.2)
  summ <- rbind(summ, data.frame(comparison = cp$name, n_treated = length(g_num), n_control = length(g_den), genes_tested = nrow(pl),
                  sig_padj = sum(pl$padj < alpha), up_padj = sum(pl$padj < alpha & pl$log2FoldChange > 0), down_padj = sum(pl$padj < alpha & pl$log2FoldChange < 0),
                  up_padj_lfc = tn[["Up"]], down_padj_lfc = tn[["Down"]]))
  cat(cp$name, ": tested", nrow(pl), "; padj<", alpha, ":", sum(pl$padj < alpha), "; up/down with |log2FC|>=", lfc_cut, ":", tn[["Up"]], "/", tn[["Down"]], "\n")
}
write.csv(summ, file.path(out_dir, "tables", "summary_counts.csv"), row.names = FALSE)
sheets <- c(list(summary = summ), sheets)
write.xlsx(sheets, file.path(out_dir, "tables", "DE_all_comparisons.xlsx"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nDONE\n")
