# Differential m6A (DART-seq) per site and per gene, treated vs control group, following the Meyer-lab method (Flamand 2022 Genes Dev; Tegowski 2024 Nat Neurosci).
# Usage: Rscript 04_diff_m6a.R <sites_dir> <comparisons.tsv> <out_dir> [<diff_expr_tables_dir> | -]
#   sites_dir       : output of 03_matrix.py (<cmp>/sites_kept.tsv, <cmp>/counts.tsv)
#   comparisons.tsv : columns name, treated, control [, de_name]; DE table = <diff_expr_tables_dir>/<de_name or name>.csv (deseq2-expression skill);
#                     '-' or a missing table = no expression columns
# Internal variable names: trt = treated group, ctl = control group.
# Two versions per comparison: "norac" (no RAC filter) and "rac" (only sites with R-A-C, in_rac == 1).
# Per site : keep samples with cov >= MIN_COV; skip sites with < MIN_PER_GROUP such samples in either group;
#            fit glm(mut/cov ~ condition, weights = cov, family = quasibinomial), condition reference = control, so the coefficient is the treatment effect.
#            P = Wald test of `condition` from summary(glm): for a quasi-binomial fit R reports a t statistic (estimate/SE, SE scaled by the estimated dispersion) referred to t(df.residual) = `wald_p`.
#            The paper (Tegowski 2024 Nat Neurosci) says "Wald's test" without stating t or z; `wald_z_p` (same statistic vs the normal) is given for comparison only. IHW uses wald_p.
#            Adjusted P: IHW with mean coverage as the covariate (OUR choice, not from the papers); fallback to BH only if IHW fails, and the column says which was used.
#            IHW splits the tests into random folds: set.seed(SEED) before every ihw() call makes the adjusted P reproducible.
# Per gene : sum mut and cov over all kept sites of the gene per sample, same GLM; RNA editing score = sum(mut)/sum(cov) * n_sites; genes with >= GENE_FC-fold change and adj. P < ALPHA.
# Known limitations (kept on purpose, see SKILL.md): no random effect for replicate (overdispersion only through the quasi-binomial dispersion);
#   when one group has 0 edited reads at a site the Wald test is uninformative (Hauck-Donner, P ~ 1) although the difference can be large.
# env (optional): MIN_COV (20), MIN_PER_GROUP (2), ALPHA (0.05), GENE_FC (1.25), EPS (0.01, pseudo-frequency for the site log2FC), SEED (1)
suppressPackageStartupMessages({library(data.table); library(IHW); library(ggplot2); library(ggrepel)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript 04_diff_m6a.R <sites_dir> <comparisons.tsv> <out_dir> [<diff_expr_tables_dir> | -]")
sites_dir <- args[1]; cmp_file <- args[2]; out_dir <- args[3]; de_dir <- if (length(args) >= 4) args[4] else "-"
envn <- function(k, d) as.numeric(Sys.getenv(k, d))
MIN_COV <- envn("MIN_COV", "20"); MIN_PER_GROUP <- envn("MIN_PER_GROUP", "2"); ALPHA <- envn("ALPHA", "0.05"); GENE_FC <- envn("GENE_FC", "1.25")
EPS <- envn("EPS", "0.01"); SEED <- as.integer(envn("SEED", "1"))
cat("MIN_COV", MIN_COV, " MIN_PER_GROUP", MIN_PER_GROUP, " ALPHA", ALPHA, " GENE_FC", GENE_FC, " EPS", EPS, " SEED", SEED, "\n")
cf <- fread(cmp_file, colClasses = "character")
CMPS <- data.table(cmp = cf$name, trt = cf$treated, ctl = cf$control,
                   de = if ("de_name" %in% names(cf)) ifelse(is.na(cf$de_name) | cf$de_name == "", cf$name, cf$de_name) else cf$name)

fit_one <- function(mut, cov, grp) {            # grp: factor, level 1 = control, level 2 = treated; returns c(beta, p, pz) or NA
  df <- data.frame(mut = mut, cov = cov, g = grp)
  fit <- tryCatch(suppressWarnings(glm(mut / cov ~ g, data = df, weights = cov, family = quasibinomial)), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, 3))
  co <- summary(fit)$coefficients
  if (nrow(co) < 2) return(rep(NA_real_, 3))
  c(co[2, 1], co[2, 4], 2 * pnorm(-abs(co[2, 1] / co[2, 2])))   # coef, P from the t test reported by summary(), P if the same statistic were referred to the normal (reported only for comparison)
}
adjust_ihw <- function(p, cov) {                 # returns list(adj, method)
  ok <- !is.na(p)
  adj <- rep(NA_real_, length(p)); meth <- "IHW"
  set.seed(SEED)
  r <- tryCatch(ihw(p[ok], cov[ok], alpha = ALPHA), error = function(e) NULL)
  if (is.null(r)) { adj[ok] <- p.adjust(p[ok], "BH"); meth <- "BH (IHW failed)" } else adj[ok] <- adj_pvalues(r)
  list(adj = adj, method = meth)
}
test_table <- function(d, key_cols, trt, ctl) {  # d: long table with key_cols, condition, mut, cov (samples with cov>=MIN_COV only)
  d[, g := factor(condition, levels = c(ctl, trt))]
  setkeyv(d, key_cols)
  res <- d[, {
    n1 <- sum(g == ctl); n2 <- sum(g == trt)
    if (n1 < MIN_PER_GROUP || n2 < MIN_PER_GROUP) list(n_ctl = n1, n_trt = n2, tested = FALSE, beta = NA_real_, p = NA_real_, pz = NA_real_,
         f_ctl = NA_real_, f_trt = NA_real_, mean_cov = NA_real_)
    else { ft <- fit_one(mut, cov, g)
      list(n_ctl = n1, n_trt = n2, tested = TRUE, beta = ft[1], p = ft[2], pz = ft[3],
           f_ctl = sum(mut[g == ctl]) / sum(cov[g == ctl]), f_trt = sum(mut[g == trt]) / sum(cov[g == trt]), mean_cov = mean(cov)) }
  }, by = key_cols]
  res
}
add_effect <- function(res, trt, ctl) {
  res[, `:=`(diff_edit = f_trt - f_ctl, log2FC = log2((f_trt + EPS) / (f_ctl + EPS)))]
  res
}
call_dir <- function(adj, diff) ifelse(!is.na(adj) & adj < ALPHA & diff < 0, "lost_in_treated",
                                 ifelse(!is.na(adj) & adj < ALPHA & diff > 0, "gained_in_treated", "ns"))
volcano <- function(res, xcol, xlab, title, file, label_n = 15) {
  pl <- res[!is.na(adj_p)]
  pl[, y := -log10(pmax(adj_p, 1e-300))]
  pl[, dir := factor(direction, levels = c("lost_in_treated", "gained_in_treated", "ns"))]
  lab <- pl[direction != "ns"][order(adj_p)][seq_len(min(label_n, .N))]
  gp <- ggplot(pl, aes(x = .data[[xcol]], y = y, color = dir)) + geom_point(size = 0.8, alpha = 0.6) +
    scale_color_manual(values = c(lost_in_treated = "#2c7fb8", gained_in_treated = "#d95f0e", ns = "grey70"), drop = FALSE, name = NULL) +
    geom_hline(yintercept = -log10(ALPHA), linetype = 2, color = "grey40") + geom_vline(xintercept = 0, color = "grey60") +
    geom_text_repel(data = lab, aes(label = gene), size = 2.6, max.overlaps = 30, show.legend = FALSE, color = "black") +
    labs(x = xlab, y = expression(-log[10]~"adj. P"), title = title) + theme_bw(base_size = 11)
  ggsave(paste0(file, ".png"), gp, width = 6.5, height = 5, dpi = 200); ggsave(paste0(file, ".pdf"), gp, width = 6.5, height = 5)
}

summ <- list(); gsumm <- list(); coin <- list()
for (i in seq_len(nrow(CMPS))) {
  cm <- CMPS$cmp[i]; trt <- CMPS$trt[i]; ctl <- CMPS$ctl[i]
  sk <- fread(file.path(sites_dir, cm, "sites_kept.tsv"), colClasses = c(chr = "character"))
  ct <- fread(file.path(sites_dir, cm, "counts.tsv"), colClasses = c(chr = "character"))
  sk[, site := paste0(chr, ":", pos, ":", strand)]; ct[, site := paste0(chr, ":", pos, ":", strand)]
  de_file <- file.path(de_dir, paste0(CMPS$de[i], ".csv"))
  de <- if (de_dir != "-" && file.exists(de_file)) fread(de_file) else data.table(symbol = character(), log2FoldChange = numeric(), padj = numeric())
  if (de_dir != "-" && !file.exists(de_file)) cat("note: no DE table", de_file, "- expression columns left empty\n")
  for (ver in c("norac", "rac")) {
    od <- file.path(out_dir, ver, cm); dir.create(od, recursive = TRUE, showWarnings = FALSE)
    s <- if (ver == "rac") sk[in_rac == 1] else sk
    cc <- ct[site %in% s$site]
    n_sites <- nrow(s)
    wide_cov <- dcast(cc, site ~ sample, value.var = "cov")
    # ---- per-site test
    use <- cc[cov >= MIN_COV]
    res <- test_table(copy(use[, .(site, condition, mut, cov)]), "site", trt, ctl)
    res <- merge(s[, .(site, chr, pos, strand, gene)], res, by = "site", all.x = TRUE)
    res[is.na(tested), `:=`(tested = FALSE, n_ctl = 0L, n_trt = 0L)]
    n_dropped <- sum(!res$tested)
    res <- add_effect(res, trt, ctl)
    res[, adj_p := NA_real_]
    a <- adjust_ihw(res$p, res$mean_cov); res[, adj_p := a$adj]; adj_method <- a$method
    setnames(res, c("f_ctl", "f_trt", "p", "pz", "beta", "n_ctl", "n_trt"), c(paste0("edit_", ctl), paste0("edit_", trt), "wald_p", "wald_z_p", "glm_coef_logit", paste0("n_samples_", ctl), paste0("n_samples_", trt)))
    res[, direction := call_dir(adj_p, diff_edit)]
    res[, adj_method := adj_method]
    # per-sample coverage and edit frequency of the kept site, for the table
    cc[, f := ifelse(cov > 0, mut / cov, NA_real_)]
    wf <- dcast(cc, site ~ sample, value.var = c("mut", "cov"))
    out <- merge(res, wf, by = "site", all.x = TRUE)
    setcolorder(out, c("site", "gene", "chr", "pos", "strand"))
    out <- out[order(adj_p, wald_p)]
    fwrite(out, file.path(od, "per_site_results.tsv"), sep = "\t")
    volcano(res[tested == TRUE], "diff_edit", paste0("Editing frequency difference (", trt, " - ", ctl, ")"),
            paste0(cm, " [", ver, "]: ", sum(res$direction == "lost_in_treated"), " lost, ", sum(res$direction == "gained_in_treated"), " gained"), file.path(od, "volcano_sites"))
    # ---- gene level
    sg <- merge(cc, s[, .(site, gene)], by = "site")
    gs <- sg[, .(mut = sum(mut), cov = sum(cov), n_sites = uniqueN(site)), by = .(gene, sample, condition)]
    gs <- gs[cov >= MIN_COV]
    nsite_gene <- s[, .(n_kept_sites = .N), by = gene]
    gres <- test_table(copy(gs[, .(gene, condition, mut, cov)]), "gene", trt, ctl)
    gres <- merge(gres, nsite_gene, by = "gene")
    gres[, `:=`(score_ctl = f_ctl * n_kept_sites, score_trt = f_trt * n_kept_sites)]   # RNA editing score = sum(mut)/sum(cov) x number of kept sites
    gres[, `:=`(fold_change = score_trt / score_ctl, diff_edit = f_trt - f_ctl)]
    gres[, log2FC := log2(fold_change)]
    ga <- adjust_ihw(gres$p, gres$mean_cov); gres[, adj_p := ga$adj]
    gres[, direction := ifelse(!is.na(adj_p) & adj_p < ALPHA & !is.na(fold_change) & fold_change >= GENE_FC, "gained_in_treated",
                        ifelse(!is.na(adj_p) & adj_p < ALPHA & !is.na(fold_change) & fold_change <= 1 / GENE_FC, "lost_in_treated", "ns"))]
    gres[, adj_method := ga$method]
    setnames(gres, c("p", "pz", "beta"), c("wald_p", "wald_z_p", "glm_coef_logit"))
    # ---- check against RNA expression (deseq2-expression table of the same comparison, control-sample counts)
    gres <- merge(gres, de[, .(gene = symbol, expr_log2FC = log2FoldChange, expr_padj = padj)][!duplicated(gene)], by = "gene", all.x = TRUE)
    gres <- gres[order(adj_p, wald_p)]
    fwrite(gres, file.path(od, "per_gene_results.tsv"), sep = "\t")
    gv <- copy(gres)
    volcano(gv[tested == TRUE], "log2FC", paste0("log2 fold change of RNA editing score (", trt, " / ", ctl, ")"),
            paste0(cm, " [", ver, "] gene level"), file.path(od, "volcano_genes"))
    # expression coincidence of differential sites (by gene)
    sig <- out[direction != "ns"]
    sg_gene <- unique(sig[, .(gene, direction)])
    mm <- merge(sig[, .(gene)][, .N, by = gene], de[, .(gene = symbol, expr_log2FC = log2FoldChange, expr_padj = padj)][!duplicated(gene)], by = "gene", all.x = TRUE)
    n_sig_genes <- nrow(mm); n_expr_sig <- sum(!is.na(mm$expr_padj) & mm$expr_padj < ALPHA)
    summ[[length(summ) + 1]] <- data.table(comparison = cm, version = ver, sites_kept = n_sites, sites_tested = sum(res$tested),
      sites_dropped_cov_filter = n_dropped, sites_lost_in_treated = sum(res$direction == "lost_in_treated"), sites_gained_in_treated = sum(res$direction == "gained_in_treated"),
      adj_method = adj_method, genes_tested = sum(gres$tested), genes_lost = sum(gres$direction == "lost_in_treated"), genes_gained = sum(gres$direction == "gained_in_treated"),
      genes_with_diff_sites = n_sig_genes, of_which_expression_padj_lt_0.05 = n_expr_sig,
      pct_expr_sig = round(100 * n_expr_sig / max(n_sig_genes, 1), 1),
      background_pct_expr_sig = round(100 * mean(de$padj < ALPHA, na.rm = TRUE), 1))
    cat(cm, ver, ": kept", n_sites, " tested", sum(res$tested), " dropped(cov)", n_dropped, " lost", sum(res$direction == "lost_in_treated"),
        " gained", sum(res$direction == "gained_in_treated"), " [", adj_method, "]\n", sep = "")
  }
}
summ <- rbindlist(summ)
fwrite(summ, file.path(out_dir, "summary_counts.tsv"), sep = "\t")
print(summ)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
