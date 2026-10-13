#!/usr/bin/env Rscript
# Extra numbers and two figures for the dataset-2 slides (Indian Liver Patient Records).
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- if (length(this_file)) normalizePath(file.path(dirname(this_file[1]), "..")) else getwd()
source(file.path(root, "scripts", "utils", "theme.R"))
out <- file.path(root, "presentation", "figures"); res <- file.path(root, "presentation", "results")
skew <- function(x) { m <- mean(x); mean((x - m)^3) / (mean((x - m)^2))^1.5 }
savep <- function(p, name, w, h) ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, dpi = 220, bg = "white",
                                        device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else grDevices::png)
liver <- read.csv(file.path(root, "data/raw/indian_liver_patient.csv")) |>
  mutate(status = factor(Dataset, levels = c(1, 2), labels = c("Liver patient", "No liver disease")))

# ---- normality + skewness for age and the 8 biomarkers (all 583 patients) ---------------
vars <- c(Age = "Age (years)", Total_Bilirubin = "Total bilirubin", Direct_Bilirubin = "Direct bilirubin",
          Alkaline_Phosphotase = "Alkaline phosphatase", Alamine_Aminotransferase = "ALT",
          Aspartate_Aminotransferase = "AST", Total_Protiens = "Total protein", Albumin = "Albumin",
          Albumin_and_Globulin_Ratio = "Albumin/globulin ratio")
norm_tab <- lapply(names(vars), function(v) { x <- na.omit(liver[[v]]); s <- shapiro.test(x)
  data.frame(variable = vars[[v]], n = length(x), mean = mean(x), median = median(x), skewness = skew(x),
             W = s$statistic, p = s$p.value) }) |> bind_rows()
write.csv(norm_tab, file.path(res, "liver_normality_skewness.csv"), row.names = FALSE); print(norm_tab, digits = 3)

# ---- box plots of five tests by diagnosis, Wilcoxon p in panel titles ---------------------
tests <- c(Total_Bilirubin = "Total bilirubin (mg/dL, log10)", Alkaline_Phosphotase = "Alkaline phosphatase (IU/L, log10)",
           Alamine_Aminotransferase = "ALT (IU/L, log10)", Aspartate_Aminotransferase = "AST (IU/L, log10)", Albumin = "Albumin (g/dL)")
b <- liver |> transmute(status, Total_Bilirubin = log10(Total_Bilirubin), Alkaline_Phosphotase = log10(Alkaline_Phosphotase),
                        Alamine_Aminotransferase = log10(Alamine_Aminotransferase), Aspartate_Aminotransferase = log10(Aspartate_Aminotransferase),
                        Albumin) |>
  pivot_longer(-status, names_to = "test") |>
  mutate(test = factor(tests[test], levels = tests))
wp <- lapply(names(tests), function(v) data.frame(test = tests[[v]], p = wilcox.test(liver[[v]] ~ liver$status)$p.value)) |> bind_rows() |>
  mutate(test = factor(test, levels = tests), label = paste("Wilcoxon", fmt_p(p)))
n_st <- liver |> count(status)
fig_box <- ggplot(b, aes(status, value, fill = status)) +
  geom_boxplot(width = 0.55, alpha = 0.6, outlier.shape = 21, outlier.size = 1.2) +
  geom_text(data = wp, aes(x = 0.5, y = Inf, label = label), hjust = 0, vjust = 1.5, size = 3.4, colour = "#52514e", inherit.aes = FALSE) +
  facet_wrap(~ test, scales = "free_y", nrow = 1, labeller = labeller(test = function(x) sub(" \\(", "\n(", x))) +
  scale_fill_manual(values = pal_liver, guide = "none") +
  scale_x_discrete(labels = setNames(sprintf("%s\n(n = %d)", c("Liver\npatients", "No liver\ndisease"), n_st$n), n_st$status)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(x = NULL, y = NULL, title = "Liver function tests are higher, and albumin lower, in liver patients",
       subtitle = "Box = median and IQR, whiskers = 1.5 x IQR, hollow points = outliers; enzymes and bilirubin on a log10 scale. Two-sided Wilcoxon rank-sum tests") +
  theme(panel.grid.major.x = element_blank())
savep(fig_box, "liver_fig_boxplots_wilcoxon", 12.5, 5.6)

# ---- regression diagnostics for Albumin ~ Total protein -------------------------------------
fit <- lm(Albumin ~ Total_Protiens, data = liver)
d <- data.frame(fitted = fitted(fit), resid = resid(fit), id = seq_len(nrow(liver)), status = liver$status)
big <- d |> arrange(resid) |> head(3)
p1 <- ggplot(d, aes(fitted, resid)) +
  geom_hline(yintercept = 0, colour = grey_mark) +
  geom_point(aes(colour = status), size = 1.6, alpha = 0.7) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.8) +
  geom_text(data = big, aes(label = paste("patient", id)), vjust = -0.8, size = 3.2, colour = "#52514e") +
  scale_colour_manual(values = pal_liver, name = NULL) +
  labs(x = "Fitted albumin (g/dL)", y = "Residual (observed - fitted, g/dL)",
       title = "Residuals vs fitted values", subtitle = "Black curve = smoothed trend; flat around zero = linear model adequate")
sw <- shapiro.test(d$resid)
p2 <- ggplot(d, aes(sample = resid)) +
  geom_abline(intercept = mean(d$resid), slope = sd(d$resid), colour = grey_mark) +
  stat_qq(colour = pal_liver[["Liver patient"]], size = 1.4, alpha = 0.7) +
  labs(x = "Theoretical quantiles (standard normal)", y = "Residual quantiles",
       title = "Normal Q-Q plot of the residuals",
       subtitle = sprintf("Shapiro-Wilk W = %.3f, %s; a few large negative residuals", sw$statistic, fmt_p(sw$p.value)))
library(patchwork)
fig_diag <- (p1 | p2) + plot_annotation(title = "Checking the assumptions of the regression Albumin ~ Total protein (n = 583)",
                                        theme = theme(plot.title = element_text(face = "bold", size = 15)))
savep(fig_diag, "liver_fig_regression_diagnostics", 13, 5.6)
writeLines(c(capture.output(print(summary(fit))), sprintf("Shapiro-Wilk on residuals: W = %.4f, p = %.3g", sw$statistic, sw$p.value),
             capture.output(print(liver[big$id, c("Total_Protiens", "Albumin", "Albumin_and_Globulin_Ratio", "Dataset")]))),
           file.path(res, "liver_regression_diagnostics.txt"))
message("liver extras done")
