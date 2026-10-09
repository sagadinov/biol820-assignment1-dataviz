#!/usr/bin/env Rscript
# ===========================================================================
# Dataset 1 - Breast cancer gene expression (CuMiDa, GEO series GSE45827)
# BIOL 520/820 - Assignment 1 (Data Visualization)
#
# Kaggle: https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida
# 151 samples x 54,675 Affymetrix HG-U133 Plus 2.0 probes, 6 classes
# (normal, luminal A, luminal B, HER2-enriched, basal-like, cell line).
#
# Figures produced (figures/pdf + figures/png):
#   fig1_1_pca            PCA score plot + scree bars        (scatter-plot family)
#   fig1_2_heatmap        clustered heatmap, top-500 variable probes (heatmap)
#   fig1_3_volcano        volcano plot, basal-like vs luminal A (scatter-plot family)
#   fig1_4_marker_violins violin + box plots of marker genes by subtype (box/violin)
# Tables produced (results/): class counts, PCA variance, DE results, marker probes.
# ===========================================================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(ggplot2)
  library(ggrepel); library(patchwork); library(pheatmap); library(matrixStats)
})
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- if (length(this_file)) normalizePath(file.path(dirname(this_file[1]), "..")) else normalizePath(getwd())
source(file.path(root, "scripts", "utils", "theme.R"))
res_dir <- file.path(root, "results"); dir.create(res_dir, showWarnings = FALSE)
set.seed(820)

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
message("Reading expression matrix ...")
dt <- fread(file.path(root, "data/raw/Breast_GSE45827.csv"))
stopifnot(all(c("samples", "type") %in% names(dt)))
expr <- t(as.matrix(dt[, -(1:2)]))          # probes x samples (log2 RMA values)
colnames(expr) <- dt$samples
subtype <- factor(dt$type, levels = subtype_levels, labels = subtype_labels)
message(sprintf("%d probes x %d samples; expression range %.2f - %.2f",
                nrow(expr), ncol(expr), min(expr), max(expr)))

class_counts <- data.frame(subtype = levels(subtype), n = as.integer(table(subtype)))
fwrite(class_counts, file.path(res_dir, "ds1_class_counts.csv"))
fwrite(data.frame(metric = c("n_samples", "n_probes", "min_expr", "max_expr", "median_expr"),
                  value = c(ncol(expr), nrow(expr), min(expr), max(expr), median(expr))),
       file.path(res_dir, "ds1_dataset_info.csv"))

# ---- probe -> gene symbol annotation (GPL570, NCBI GEO) --------------------
ann_lines <- readLines(gzfile(file.path(root, "data/raw/GPL570.annot.gz")))
b <- grep("^!platform_table_begin", ann_lines); e <- grep("^!platform_table_end", ann_lines)
ann <- fread(text = ann_lines[(b + 1):(e - 1)], sep = "\t", quote = "",
             select = c("ID", "Gene symbol", "Gene title"))
setnames(ann, c("probe", "symbol", "title"))
sym_full  <- setNames(ann$symbol, ann$probe)
sym_first <- setNames(sub("///.*$", "", ann$symbol), ann$probe)   # first symbol if several
symbol_of <- function(p) { s <- sym_first[p]; s[is.na(s)] <- ""; unname(s) }

# Pick, for a gene symbol, the probe with the highest mean expression.
probe_for <- function(gene) {
  hits <- ann$probe[vapply(strsplit(ann$symbol, "///", fixed = TRUE),
                           function(s) gene %in% s, logical(1))]
  hits <- intersect(hits, rownames(expr))
  if (!length(hits)) return(NA_character_)
  hits[which.max(rowMeans(expr[hits, , drop = FALSE]))]
}

# ---------------------------------------------------------------------------
# 2. Figure 1.1 - PCA score plot (top 2,000 most variable probes) + scree
# ---------------------------------------------------------------------------
rv  <- rowVars(expr)
top2000 <- order(rv, decreasing = TRUE)[1:2000]
pca <- prcomp(t(expr[top2000, ]), center = TRUE, scale. = FALSE)
ve  <- pca$sdev^2 / sum(pca$sdev^2)
fwrite(data.frame(PC = paste0("PC", seq_along(ve))[1:10], variance_explained = ve[1:10],
                  cumulative = cumsum(ve)[1:10]), file.path(res_dir, "ds1_pca_variance.csv"))

scores <- data.frame(pca$x[, 1:4], sample = colnames(expr), subtype = subtype)
fwrite(scores, file.path(res_dir, "ds1_pca_scores.csv"))
centroids <- scores |> group_by(subtype) |>
  summarise(PC1 = median(PC1), PC2 = median(PC2), n = n(), .groups = "drop") |>
  mutate(label = sprintf("%s (n = %d)", subtype, n))

p_pca <- ggplot(scores, aes(PC1, PC2)) +
  stat_ellipse(aes(colour = subtype), type = "norm", level = 0.95, linewidth = 0.4, alpha = 0.9,
               show.legend = FALSE) +
  geom_point(colour = ink$surface, size = 3.4) +                      # 2-px surface ring
  geom_point(aes(colour = subtype, shape = subtype), size = 2.4) +
  geom_label_repel(data = centroids, aes(label = label), size = 3, colour = ink$primary,
                   fill = scales::alpha(ink$surface, 0.85), label.size = 0, seed = 1,
                   min.segment.length = 0, segment.colour = ink$muted, box.padding = 0.6) +
  scale_colour_manual(values = pal_subtype, name = NULL) +
  scale_shape_manual(values = shape_subtype, name = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3))) +
  labs(x = sprintf("PC1 (%.1f%% of variance)", 100 * ve[1]),
       y = sprintf("PC2 (%.1f%% of variance)", 100 * ve[2]),
       title = "Breast tumour samples separate by molecular subtype",
       subtitle = "PCA of the 2,000 most variable probes, GSE45827 (n = 151); ellipses = 95% normal-probability contours") +
  theme(legend.position = "top")

scree <- data.frame(PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
                    ve = 100 * ve[1:10], shown = c(TRUE, TRUE, rep(FALSE, 8)))
p_scree <- ggplot(scree, aes(PC, ve, fill = shown)) +
  geom_col(width = 0.65) +
  geom_text(data = subset(scree, shown), aes(label = sprintf("%.1f%%", ve)), vjust = -0.4,
            size = 2.8, colour = ink$primary) +
  scale_fill_manual(values = c(`TRUE` = pal_cat[["blue"]], `FALSE` = col_ns), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Variance explained (%)", title = "Scree plot") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = rel(1)))

fig1_1 <- p_pca + p_scree + plot_layout(widths = c(3.2, 1)) +
  plot_annotation(caption = "Data: CuMiDa Breast_GSE45827 (Feltes et al., 2019; Gruosso et al., 2016). Log2 RMA-normalised Affymetrix HG-U133 Plus 2.0 values.",
                  theme = theme(plot.caption = element_text(colour = ink$muted, size = 8, hjust = 0)))
save_fig(fig1_1, "fig1_1_pca", root, width = 11.5, height = 6.2)

# ---------------------------------------------------------------------------
# 3. Figure 1.2 - Clustered heatmap of the 500 most variable probes
# ---------------------------------------------------------------------------
top500 <- order(rv, decreasing = TRUE)[1:500]
z <- t(scale(t(expr[top500, ])))                  # row z-scores
z[z >  3] <-  3; z[z < -3] <- -3

markers_hm <- c("ESR1", "PGR", "FOXA1", "GATA3", "XBP1", "TFF1", "AGR2", "SCGB2A2",
                "ERBB2", "KRT5", "KRT14", "KRT17", "FOXC1", "EGFR", "SFRP1",
                "MKI67", "TOP2A", "BIRC5", "CDH1", "VIM", "KIT", "ELF5", "MIA", "GABRP",
                "S100A8", "CA12", "NAT1", "SLC39A6", "ANKRD30A", "CDC20", "CCNB1")
row_lab <- symbol_of(rownames(z)); row_lab[!(row_lab %in% markers_hm)] <- ""
row_lab[duplicated(row_lab) & row_lab != ""] <- ""     # label each marker once

esr1_p  <- probe_for("ESR1"); erbb2_p <- probe_for("ERBB2"); krt5_p <- probe_for("KRT5")
ann_col <- data.frame(ESR1 = expr[esr1_p, ], ERBB2 = expr[erbb2_p, ], KRT5 = expr[krt5_p, ],
                      Subtype = subtype, row.names = colnames(expr))
seq_ramp <- colorRampPalette(pal_seq)(50)
ann_colors <- list(Subtype = pal_subtype, ESR1 = seq_ramp, ERBB2 = seq_ramp, KRT5 = seq_ramp)
col_lab <- c("Normal" = "N", "Luminal A" = "LA", "Luminal B" = "LB", "HER2-enriched" = "H2",
             "Basal-like" = "B", "Cell line" = "CL")[as.character(subtype)]

ph <- pheatmap(z,
               color = colorRampPalette(c(pal_div[["low"]], pal_div[["mid"]], pal_div[["high"]]))(101),
               breaks = seq(-3, 3, length.out = 102),
               clustering_distance_cols = "correlation", clustering_distance_rows = "euclidean",
               clustering_method = "ward.D2",
               annotation_col = ann_col, annotation_colors = ann_colors,
               labels_row = row_lab, labels_col = col_lab,
               fontsize = 8, fontsize_row = 6.5, fontsize_col = 4.5,
               treeheight_row = 22, treeheight_col = 40, border_color = NA,
               main = "Unsupervised clustering of the 500 most variable probes recovers the breast cancer subtypes (row z-scores)",
               silent = TRUE)
save_grob(ph$gtable, "fig1_2_heatmap", root, width = 12, height = 9)

# Which subtype dominates each of k = 6 column clusters? (for the discussion)
col_clusters <- cutree(ph$tree_col, k = 6)
cluster_tab <- as.data.frame.matrix(table(cluster = col_clusters, subtype = subtype))
cluster_tab <- cbind(cluster = rownames(cluster_tab), cluster_tab)
fwrite(cluster_tab, file.path(res_dir, "ds1_heatmap_column_clusters.csv"))

# ---------------------------------------------------------------------------
# 4. Figure 1.3 - Volcano plot: basal-like vs luminal A (Welch t-test per probe)
# ---------------------------------------------------------------------------
g1 <- subtype == "Basal-like"; g2 <- subtype == "Luminal A"
m1 <- rowMeans(expr[, g1]); m2 <- rowMeans(expr[, g2])
v1 <- rowVars(expr[, g1]);  v2 <- rowVars(expr[, g2])
n1 <- sum(g1); n2 <- sum(g2)
se2 <- v1 / n1 + v2 / n2
tstat <- (m1 - m2) / sqrt(se2)
df_w  <- se2^2 / ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))
pval  <- 2 * pt(-abs(tstat), df_w)
de <- data.frame(probe = rownames(expr), symbol = symbol_of(rownames(expr)),
                 mean_basal = m1, mean_luminalA = m2, log2FC = m1 - m2,
                 t = tstat, df = df_w, p = pval, stringsAsFactors = FALSE) |>
  filter(is.finite(p)) |>
  mutate(fdr = p.adjust(p, method = "BH"),
         status = case_when(fdr < 0.05 & log2FC >=  1 ~ "Higher in basal-like",
                            fdr < 0.05 & log2FC <= -1 ~ "Higher in luminal A",
                            TRUE ~ "Not significant"),
         status = factor(status, levels = c("Higher in basal-like", "Higher in luminal A", "Not significant")))
n_up <- sum(de$status == "Higher in basal-like"); n_down <- sum(de$status == "Higher in luminal A")
p_cut <- max(de$p[de$fdr < 0.05])                   # raw p corresponding to FDR = 0.05
fwrite(de |> filter(status != "Not significant") |> arrange(fdr),
       file.path(res_dir, "ds1_de_basal_vs_luminalA_significant.csv"))
fwrite(data.frame(metric = c("n_tested", "n_higher_basal", "n_higher_luminalA", "p_at_fdr05", "n_basal", "n_luminalA"),
                  value = c(nrow(de), n_up, n_down, p_cut, n1, n2)),
       file.path(res_dir, "ds1_de_summary.csv"))

lab_df <- de |> filter(status != "Not significant", symbol != "") |>
  group_by(status) |> arrange(fdr, .by_group = TRUE) |>
  distinct(symbol, .keep_all = TRUE) |> slice_head(n = 12) |> ungroup()

fig1_3 <- ggplot(de, aes(log2FC, -log10(p))) +
  geom_point(data = subset(de, status == "Not significant"), colour = col_ns, size = 0.6, alpha = 0.5) +
  geom_point(data = subset(de, status != "Not significant"), aes(colour = status), size = 0.9, alpha = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "22", colour = ink$muted, linewidth = 0.35) +
  geom_hline(yintercept = -log10(p_cut), linetype = "22", colour = ink$muted, linewidth = 0.35) +
  geom_text_repel(data = lab_df, aes(label = symbol), size = 2.8, colour = ink$primary,
                  max.overlaps = Inf, box.padding = 0.35, min.segment.length = 0,
                  segment.colour = ink$muted, segment.size = 0.25, seed = 2) +
  annotate("label", x = max(de$log2FC), y = 0.5, hjust = 1, vjust = 0, linewidth = 0,
           label = sprintf("%s probes higher in basal-like", format(n_up, big.mark = ",")),
           size = 3, colour = ink$primary, fill = scales::alpha(ink$surface, 0.8)) +
  annotate("label", x = min(de$log2FC), y = 0.5, hjust = 0, vjust = 0, linewidth = 0,
           label = sprintf("%s probes higher in luminal A", format(n_down, big.mark = ",")),
           size = 3, colour = ink$primary, fill = scales::alpha(ink$surface, 0.8)) +
  scale_colour_manual(values = c("Higher in basal-like" = pal_div[["high"]],
                                 "Higher in luminal A"  = pal_div[["low"]]), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(x = expression(log[2]~"fold change (basal-like − luminal A)"),
       y = expression(-log[10]~"(p-value)"),
       title = "Basal-like and luminal A tumours differ in thousands of transcripts",
       subtitle = sprintf("Welch t-test per probe (%d basal-like vs %d luminal A samples); dashed lines: |log2FC| = 1 and FDR = 0.05 (Benjamini–Hochberg)", n1, n2),
       caption = sprintf("%s probes tested; labelled genes are the 12 most significant in each direction. Gene symbols from the NCBI GEO GPL570 annotation.",
                         format(nrow(de), big.mark = ",")))
save_fig(fig1_3, "fig1_3_volcano", root, width = 10, height = 7)

# ---------------------------------------------------------------------------
# 5. Figure 1.4 - Marker-gene expression by subtype (violin + box + points)
# ---------------------------------------------------------------------------
marker_genes <- c("ESR1", "PGR", "ERBB2", "KRT5", "KRT17", "MKI67")
marker_probe <- vapply(marker_genes, probe_for, character(1))
fwrite(data.frame(gene = marker_genes, probe = marker_probe,
                  role = c("Oestrogen receptor (luminal marker)", "Progesterone receptor (luminal marker)",
                           "HER2 receptor (HER2-enriched marker)", "Cytokeratin 5 (basal marker)",
                           "Cytokeratin 17 (basal marker)", "Ki-67 (proliferation)")),
       file.path(res_dir, "ds1_marker_probes.csv"))

long <- as.data.frame(t(expr[marker_probe, ])) |>
  setNames(marker_genes) |>
  mutate(sample = colnames(expr), subtype = subtype) |>
  pivot_longer(all_of(marker_genes), names_to = "gene", values_to = "expr") |>
  mutate(gene = factor(gene, levels = marker_genes))
short_lab <- c("Normal" = "Normal", "Luminal A" = "Lum A", "Luminal B" = "Lum B",
               "HER2-enriched" = "HER2", "Basal-like" = "Basal", "Cell line" = "Cell line")
n_txt <- paste(sprintf("%s %d", class_counts$subtype, class_counts$n), collapse = ", ")

# Kruskal-Wallis across the five tissue classes (cell lines excluded) per gene
kw <- long |> filter(subtype != "Cell line") |> group_by(gene) |>
  summarise(p = kruskal.test(expr ~ droplevels(subtype))$p.value, .groups = "drop")
fwrite(kw, file.path(res_dir, "ds1_marker_kruskal.csv"))

fig1_4 <- ggplot(long, aes(subtype, expr)) +
  geom_violin(aes(fill = subtype), colour = NA, alpha = 0.35, width = 0.9, trim = TRUE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = ink$surface, colour = ink$secondary,
               linewidth = 0.35) +
  geom_jitter(aes(colour = subtype), width = 0.14, size = 0.7, alpha = 0.6) +
  geom_text(data = kw, aes(x = 0.6, y = Inf, label = paste("Kruskal–Wallis", fmt_p(p))),
            hjust = 0, vjust = 1.5, size = 2.6, colour = ink$secondary, inherit.aes = FALSE) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = pal_subtype, guide = "none") +
  scale_colour_manual(values = pal_subtype, guide = "none") +
  scale_x_discrete(labels = short_lab) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = expression(log[2]~"expression (RMA)"),
       title = "Canonical marker genes track the clinical breast cancer subtypes",
       subtitle = sprintf("Violin = distribution, box = median and IQR, points = individual samples\nSample sizes: %s", n_txt),
       caption = sprintf("One probe per gene (highest mean expression): %s.\nKruskal–Wallis tests compare the five tissue classes (cell lines excluded).",
                         paste(sprintf("%s = %s", marker_genes, marker_probe), collapse = ", "))) +
  theme(axis.text.x = element_text(size = 8), panel.grid.major.x = element_blank())
save_fig(fig1_4, "fig1_4_marker_violins", root, width = 11.5, height = 7.2)

writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_ds1.txt"))
message("Dataset 1 done.")
