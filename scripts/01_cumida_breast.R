#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Dataset 1 - Breast cancer gene expression, CuMiDa GSE45827 (151 samples, 6 classes)
# Kaggle: https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida
#
# Figures: fig1_histogram_esr1, fig2_boxplot_subtype, fig3_scatter_esr1_gata3, fig4_qq_esr1
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(tidyr); library(ggplot2) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- if (length(this_file)) normalizePath(file.path(dirname(this_file[1]), "..")) else getwd()
source(file.path(root, "scripts", "utils", "theme.R"))
res <- file.path(root, "results"); dir.create(res, showWarnings = FALSE)

# ---- 1. Read the data and keep four well-known genes ------------------------
# Probe IDs (Affymetrix HG-U133 Plus 2.0) taken from the NCBI GEO GPL570 annotation.
genes <- c(ESR1 = "205225_at", GATA3 = "209602_s_at", ERBB2 = "216836_s_at", PGR = "228554_at")
raw <- fread(file.path(root, "data/raw/Breast_GSE45827.csv"), select = c("samples", "type", unname(genes)))
setnames(raw, genes, names(genes))
d <- raw |>
  mutate(subtype = factor(type, levels = c("normal", "luminal_A", "luminal_B", "HER", "basal", "cell_line"),
                          labels = names(pal_subtype))) |>
  select(sample = samples, subtype, ESR1, GATA3, ERBB2, PGR)
cat("Data frame:", nrow(d), "samples x", ncol(d), "columns\n"); str(d)
write.csv(d, file.path(res, "ds1_four_genes.csv"), row.names = FALSE)

# ---- 2. Descriptive statistics (week 2) --------------------------------------
desc_all <- d |> pivot_longer(ESR1:PGR, names_to = "gene") |> group_by(gene) |>
  summarise(n = n(), mean = mean(value), sd = sd(value), median = median(value),
            Q1 = quantile(value, .25), Q3 = quantile(value, .75), min = min(value), max = max(value),
            .groups = "drop")
write.csv(desc_all, file.path(res, "ds1_summary_all_samples.csv"), row.names = FALSE)
desc_sub <- d |> group_by(subtype) |>
  summarise(n = n(), ESR1_mean = mean(ESR1), ESR1_sd = sd(ESR1), ESR1_median = median(ESR1),
            ERBB2_mean = mean(ERBB2), ERBB2_sd = sd(ERBB2), ERBB2_median = median(ERBB2), .groups = "drop")
write.csv(desc_sub, file.path(res, "ds1_summary_by_subtype.csv"), row.names = FALSE)

# ---- Figure 1: histogram + density of ESR1 expression (all samples) ----------
m <- mean(d$ESR1); md <- median(d$ESR1); s <- sd(d$ESR1)
fig1 <- ggplot(d, aes(ESR1)) +
  geom_histogram(aes(y = after_stat(density), fill = subtype), bins = 30, colour = "white", linewidth = 0.3) +
  geom_density(linewidth = 0.9, colour = "black") +
  geom_vline(xintercept = m, linetype = "dashed") +
  geom_vline(xintercept = md, linetype = "dotted") +
  annotate("text", x = m + 0.15, y = Inf, vjust = 1.5, hjust = 0, size = 3.5,
           label = sprintf("mean = %.2f (dashed)\nmedian = %.2f (dotted)\nSD = %.2f", m, md, s)) +
  scale_fill_manual(values = pal_subtype, name = NULL) +
  labs(x = "ESR1 expression (log2, RMA)", y = "Density",
       title = "ESR1 (oestrogen receptor) expression is bimodal across the 151 samples",
       subtitle = "Histogram (30 bins, coloured by subtype) with kernel density curve")
save_fig(fig1, "fig1_histogram_esr1", root, width = 9, height = 5.5)

# ---- Figure 2: box plots of ESR1 and ERBB2 by subtype -------------------------
long2 <- d |> pivot_longer(c(ESR1, ERBB2), names_to = "gene") |>
  mutate(gene = factor(gene, levels = c("ESR1", "ERBB2"),
                       labels = c("ESR1 (oestrogen receptor)", "ERBB2 (HER2)")))
short <- c("Normal" = "Normal", "Luminal A" = "Lum A", "Luminal B" = "Lum B", "HER2-enriched" = "HER2",
           "Basal-like" = "Basal", "Cell line" = "Cell line")
n_lab <- d |> count(subtype) |> mutate(lab = sprintf("%s\n(n = %d)", short[as.character(subtype)], n))
fig2 <- ggplot(long2, aes(subtype, value, fill = subtype)) +
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = 21) +
  geom_jitter(width = 0.15, size = 0.9, alpha = 0.5) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = pal_subtype, guide = "none") +
  scale_x_discrete(labels = setNames(n_lab$lab, n_lab$subtype)) +
  labs(x = NULL, y = "Expression (log2, RMA)",
       title = "ESR1 and ERBB2 expression by breast cancer subtype",
       subtitle = "Box = median and inter-quartile range, whiskers = 1.5 x IQR, hollow points = outliers, dots = samples") +
  theme(axis.text.x = element_text(size = 10))
save_fig(fig2, "fig2_boxplot_subtype", root, width = 10, height = 5.5)

# ---- Figure 3: scatter plot ESR1 vs GATA3 with regression line ----------------
ct  <- cor.test(d$ESR1, d$GATA3)
fit <- lm(GATA3 ~ ESR1, data = d)
co  <- coef(fit); r2 <- summary(fit)$r.squared
writeLines(c(capture.output(print(ct)), "", capture.output(print(summary(fit)))),
           file.path(res, "ds1_correlation_regression_ESR1_GATA3.txt"))
fig3 <- ggplot(d, aes(ESR1, GATA3)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.8, fill = "grey80") +
  geom_point(aes(colour = subtype), size = 2, alpha = 0.85) +
  scale_colour_manual(values = pal_subtype, name = NULL) +
  labs(x = "ESR1 expression (log2, RMA)", y = "GATA3 expression (log2, RMA)",
       title = "GATA3 expression increases linearly with ESR1 expression",
       subtitle = sprintf("Pearson r = %.2f (%s), n = %d\nLeast-squares line: GATA3 = %.2f + %.2f x ESR1, R² = %.2f",
                          ct$estimate, fmt_p(ct$p.value), nrow(d), co[1], co[2], r2),
       caption = "Grey band = 95% confidence interval of the regression line.")
save_fig(fig3, "fig3_scatter_esr1_gata3", root, width = 9, height = 6)

# ---- Figure 4: normal Q-Q plots of ESR1 (all samples vs luminal A only) --------
qq <- bind_rows(d |> mutate(group = "All samples"),
                d |> filter(subtype == "Luminal A") |> mutate(group = "Luminal A only"))
sw <- qq |> group_by(group) |>
  summarise(n = n(), mean = mean(ESR1), sd = sd(ESR1),
            W = shapiro.test(ESR1)$statistic, p = shapiro.test(ESR1)$p.value, .groups = "drop") |>
  mutate(panel = sprintf("%s (n = %d): Shapiro-Wilk W = %.3f, %s", group, n, W, fmt_p(p)))
write.csv(sw, file.path(res, "ds1_shapiro_esr1.csv"), row.names = FALSE)
qq <- qq |> left_join(sw |> select(group, panel), by = "group")
fig4 <- ggplot(qq, aes(sample = ESR1)) +
  geom_abline(data = sw, aes(intercept = mean, slope = sd), colour = grey_mark) +   # normal with same mean and SD
  stat_qq(aes(colour = subtype), size = 1.8, alpha = 0.85) +
  facet_wrap(~ panel, scales = "free") +
  scale_colour_manual(values = pal_subtype, name = NULL) +
  labs(x = "Theoretical quantiles (standard normal)", y = "Sample quantiles of ESR1 expression",
       title = "Normal Q-Q plots of ESR1 expression",
       subtitle = "Line = quantiles of a normal distribution with the same mean and SD.\nThe step in the left panel reflects the two modes seen in Figure 1")
save_fig(fig4, "fig4_qq_esr1", root, width = 10, height = 5.2)

writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_ds1.txt"))
message("Dataset 1 done.")
