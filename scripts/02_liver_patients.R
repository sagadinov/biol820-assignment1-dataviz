#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Dataset 2 - Indian Liver Patient Records (ILPD): 583 patients, 11 columns
# Kaggle: "Indian Liver Patient Records" (same data as uciml/indian-liver-patient-records)
#
# Figures: fig5_histogram_bilirubin, fig6_boxplot_liver_tests, fig7_scatter_protein_albumin, fig8_qq_alt
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- if (length(this_file)) normalizePath(file.path(dirname(this_file[1]), "..")) else getwd()
source(file.path(root, "scripts", "utils", "theme.R"))
res <- file.path(root, "results"); dir.create(res, showWarnings = FALSE)

# ---- 1. Read and tidy -----------------------------------------------------------
liver <- read.csv(file.path(root, "data/raw/indian_liver_patient.csv")) |>
  mutate(Gender = factor(Gender),
         status = factor(Dataset, levels = c(1, 2), labels = names(pal_liver)))
cat("Data frame:", nrow(liver), "patients x", ncol(liver), "columns\n"); str(liver)
cat("Missing values per column:\n"); print(colSums(is.na(liver)))

# ---- 2. Descriptive statistics by group (week 2) ---------------------------------
vars <- c("Age", "Total_Bilirubin", "Direct_Bilirubin", "Alkaline_Phosphotase",
          "Alamine_Aminotransferase", "Aspartate_Aminotransferase", "Total_Protiens",
          "Albumin", "Albumin_and_Globulin_Ratio")
desc <- liver |> select(status, all_of(vars)) |>
  pivot_longer(-status, names_to = "variable") |>
  group_by(variable, status) |>
  summarise(n = sum(!is.na(value)), mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE), Q1 = quantile(value, .25, na.rm = TRUE),
            Q3 = quantile(value, .75, na.rm = TRUE), .groups = "drop")
write.csv(desc, file.path(res, "ds2_summary_by_status.csv"), row.names = FALSE)
write.csv(liver |> count(status, Gender), file.path(res, "ds2_counts_status_gender.csv"), row.names = FALSE)

# ---- Figure 5: histogram of total bilirubin, raw and log10 ------------------------
h <- liver |> transmute(status, raw = Total_Bilirubin, log10 = log10(Total_Bilirubin)) |>
  pivot_longer(c(raw, log10), names_to = "scale") |>
  mutate(scale = factor(scale, levels = c("raw", "log10"),
                        labels = c("Total bilirubin (mg/dL)", "log10 total bilirubin")))
sk <- h |> group_by(scale) |> summarise(mean = mean(value), median = median(value), .groups = "drop")
fig5 <- ggplot(h, aes(value)) +
  geom_histogram(bins = 35, fill = pal_liver[["Liver patient"]], colour = "white", linewidth = 0.3) +
  geom_vline(data = sk, aes(xintercept = mean), linetype = "dashed") +
  geom_vline(data = sk, aes(xintercept = median), linetype = "dotted") +
  facet_wrap(~ scale, scales = "free") +
  labs(x = NULL, y = "Number of patients",
       title = "Total bilirubin is strongly right-skewed; a log transformation reduces the skew",
       subtitle = "Dashed line = mean, dotted line = median (n = 583)")
save_fig(fig5, "fig5_histogram_bilirubin", root, width = 10, height = 4.8)

# ---- Figure 6: box plots of five liver tests by disease status --------------------
b <- liver |>
  transmute(status,
            `Total bilirubin\n(mg/dL, log10)` = log10(Total_Bilirubin),
            `Alkaline phosphatase\n(IU/L, log10)` = log10(Alkaline_Phosphotase),
            `ALT (IU/L, log10)` = log10(Alamine_Aminotransferase),
            `AST (IU/L, log10)` = log10(Aspartate_Aminotransferase),
            `Albumin (g/dL)` = Albumin) |>
  pivot_longer(-status, names_to = "test") |>
  mutate(test = factor(test, levels = unique(test)))
n_st <- liver |> count(status)
fig6 <- ggplot(b, aes(status, value, fill = status)) +
  geom_boxplot(width = 0.55, alpha = 0.6, outlier.shape = 21, outlier.size = 1.2) +
  facet_wrap(~ test, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal_liver, guide = "none") +
  scale_x_discrete(labels = setNames(sprintf("%s\n(n = %d)", c("Liver\npatients", "No liver\ndisease"), n_st$n), n_st$status)) +
  labs(x = NULL, y = NULL,
       title = "Liver function tests are higher, and albumin lower, in liver patients",
       subtitle = "Box = median and IQR, whiskers = 1.5 x IQR, hollow points = outliers; enzymes and bilirubin on a log10 scale")
save_fig(fig6, "fig6_boxplot_liver_tests", root, width = 11, height = 5.2)

# ---- Figure 7: scatter plot total protein vs albumin with regression line ---------
ct  <- cor.test(liver$Total_Protiens, liver$Albumin)
fit <- lm(Albumin ~ Total_Protiens, data = liver)
co  <- coef(fit); r2 <- summary(fit)$r.squared
writeLines(c(capture.output(print(ct)), "", capture.output(print(summary(fit)))),
           file.path(res, "ds2_correlation_regression_protein_albumin.txt"))
fig7 <- ggplot(liver, aes(Total_Protiens, Albumin)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.8, fill = "grey80") +
  geom_jitter(aes(colour = status), width = 0.03, height = 0.03, size = 1.7, alpha = 0.7) +
  scale_colour_manual(values = pal_liver, name = NULL) +
  labs(x = "Total protein (g/dL)", y = "Albumin (g/dL)",
       title = "Serum albumin rises linearly with total protein",
       subtitle = sprintf("Pearson r = %.2f (%s), n = %d\nLeast-squares line: Albumin = %.2f + %.2f x Total protein, R² = %.2f",
                          ct$estimate, fmt_p(ct$p.value), nrow(liver), co[1], co[2], r2),
       caption = "Points are slightly jittered because values are recorded to one decimal place. Grey band = 95% confidence interval of the line.")
save_fig(fig7, "fig7_scatter_protein_albumin", root, width = 9, height = 6)

# ---- Figure 8: normal Q-Q plots of ALT, raw and log10 ----------------------------
q <- liver |> transmute(raw = Alamine_Aminotransferase, log10 = log10(Alamine_Aminotransferase)) |>
  pivot_longer(c(raw, log10), names_to = "scale")
sw <- q |> group_by(scale) |>
  summarise(mean = mean(value), sd = sd(value),
            W = shapiro.test(value)$statistic, p = shapiro.test(value)$p.value, .groups = "drop") |>
  mutate(panel = sprintf("%s: Shapiro-Wilk W = %.3f, %s",
                         ifelse(scale == "raw", "ALT (IU/L)", "log10 ALT"), W, fmt_p(p)))
write.csv(sw, file.path(res, "ds2_shapiro_alt.csv"), row.names = FALSE)
q <- q |> left_join(sw |> select(scale, panel), by = "scale") |>
  mutate(panel = factor(panel, levels = sw$panel[c(which(sw$scale == "raw"), which(sw$scale == "log10"))]))
fig8 <- ggplot(q, aes(sample = value)) +
  geom_abline(data = sw, aes(intercept = mean, slope = sd), colour = grey_mark) +   # normal with same mean and SD
  stat_qq(colour = pal_liver[["Liver patient"]], size = 1.5, alpha = 0.7) +
  facet_wrap(~ panel, scales = "free") +
  labs(x = "Theoretical quantiles (standard normal)", y = "Sample quantiles",
       title = "Normal Q-Q plots of ALT (alanine aminotransferase) before and after log transformation",
       subtitle = "Line = quantiles of a normal distribution with the same mean and SD (n = 583).\nRaw ALT has a heavy right tail; log10 ALT is closer to, but still not, normal")
save_fig(fig8, "fig8_qq_alt", root, width = 11, height = 5)

writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_ds2.txt"))
message("Dataset 2 done.")
