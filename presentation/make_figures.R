#!/usr/bin/env Rscript
# Figures and statistics for the Group 5 presentation (14 Oct 2026).
# Dataset 1: Vertebral Column Data Set (UCI id 212 / Kaggle), 310 patients, 6 biomechanical features.
# Dataset 2: Indian Liver Patient Records (figures reused from scripts/02_liver_patients.R).
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- if (length(this_file)) normalizePath(file.path(dirname(this_file[1]), "..")) else getwd()
source(file.path(root, "scripts", "utils", "theme.R"))
out <- file.path(root, "presentation", "figures"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
res <- file.path(root, "presentation", "results"); dir.create(res, showWarnings = FALSE)
pal_vert <- c("Abnormal" = "#eb6834", "Normal" = "#2a78d6")
savep <- function(p, name, w, h) {
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, dpi = 220, bg = "white",
         device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else grDevices::png)
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h, bg = "white", device = grDevices::cairo_pdf)
  message("saved ", name)
}
skew <- function(x) { m <- mean(x); mean((x - m)^3) / (mean((x - m)^2))^1.5 }   # as in moments::skewness

# ---- data ---------------------------------------------------------------------------
f <- file.path(root, "data/raw/vertebral_column_2C.dat")
if (!file.exists(f)) stop("download https://archive.ics.uci.edu/static/public/212/vertebral+column.zip and copy column_2C.dat to data/raw/vertebral_column_2C.dat")
vert <- read.table(f, col.names = c("pelvic_incidence", "pelvic_tilt", "lumbar_lordosis_angle", "sacral_slope",
                                    "pelvic_radius", "degree_spondylolisthesis", "class")) |>
  mutate(class = factor(class, levels = c("AB", "NO"), labels = c("Abnormal", "Normal")))
feats <- names(vert)[1:6]
nice  <- c(pelvic_incidence = "Pelvic incidence", pelvic_tilt = "Pelvic tilt", lumbar_lordosis_angle = "Lumbar lordosis angle",
           sacral_slope = "Sacral slope", pelvic_radius = "Pelvic radius", degree_spondylolisthesis = "Degree of spondylolisthesis")
str(vert); print(table(vert$class))

# ---- Table: normality + skewness (all 310 patients) ---------------------------------
norm_tab <- lapply(feats, function(v) { s <- shapiro.test(vert[[v]])
  data.frame(variable = nice[[v]], mean = mean(vert[[v]]), median = median(vert[[v]]), sd = sd(vert[[v]]),
             skewness = skew(vert[[v]]), W = s$statistic, p = s$p.value) }) |> bind_rows()
write.csv(norm_tab, file.path(res, "vert_normality_skewness.csv"), row.names = FALSE); print(norm_tab)

# ---- Table: medians by group + Wilcoxon rank-sum per feature -------------------------
grp_tab <- lapply(feats, function(v) { w <- wilcox.test(vert[[v]] ~ vert$class)
  data.frame(variable = nice[[v]],
             median_abnormal = median(vert[[v]][vert$class == "Abnormal"]),
             median_normal   = median(vert[[v]][vert$class == "Normal"]),
             W = w$statistic, p = w$p.value) }) |> bind_rows()
write.csv(grp_tab, file.path(res, "vert_group_comparison_wilcoxon.csv"), row.names = FALSE); print(grp_tab)

# ---- Fig V1: box plots of the six features by class ----------------------------------
long <- vert |> pivot_longer(all_of(feats), names_to = "variable") |>
  mutate(variable = factor(nice[variable], levels = nice))
lab <- grp_tab |> mutate(variable = factor(variable, levels = nice), label = paste("Wilcoxon", fmt_p(p)))
n_cl <- vert |> count(class)
figV1 <- ggplot(long, aes(class, value, fill = class)) +
  geom_boxplot(width = 0.55, alpha = 0.6, outlier.shape = 21, outlier.size = 1.1) +
  geom_jitter(aes(colour = class), width = 0.15, size = 0.6, alpha = 0.35, show.legend = FALSE) +
  geom_text(data = lab, aes(x = 0.5, y = Inf, label = label), hjust = 0, vjust = 1.5, size = 3.1,
            colour = "#52514e", inherit.aes = FALSE) +
  facet_wrap(~ variable, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = pal_vert, guide = "none") + scale_colour_manual(values = pal_vert) +
  scale_x_discrete(labels = setNames(sprintf("%s\n(n = %d)", n_cl$class, n_cl$n), n_cl$class)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(x = NULL, y = "Value (degrees, except spondylolisthesis in mm)",
       title = "All six biomechanical measurements differ between abnormal and normal patients",
       subtitle = "Box = median and IQR, whiskers = 1.5 x IQR, hollow points = outliers; two-sided Wilcoxon rank-sum tests") +
  theme(panel.grid.major.x = element_blank())
savep(figV1, "vert_fig1_boxplots", 12, 6.4)

# ---- Fig V2: Q-Q plots of degree of spondylolisthesis by class -------------------------
sw <- vert |> group_by(class) |>
  summarise(n = n(), mean = mean(degree_spondylolisthesis), sd = sd(degree_spondylolisthesis),
            W = shapiro.test(degree_spondylolisthesis)$statistic,
            p = shapiro.test(degree_spondylolisthesis)$p.value, .groups = "drop") |>
  mutate(panel = sprintf("%s (n = %d): Shapiro-Wilk W = %.3f, %s", class, n, W, fmt_p(p)))
write.csv(sw, file.path(res, "vert_shapiro_spondylolisthesis_by_class.csv"), row.names = FALSE)
qq <- vert |> left_join(sw |> select(class, panel), by = "class")
figV2 <- ggplot(qq, aes(sample = degree_spondylolisthesis)) +
  geom_abline(data = sw, aes(intercept = mean, slope = sd), colour = grey_mark) +
  stat_qq(aes(colour = class), size = 1.6, alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~ panel, scales = "free") +
  scale_colour_manual(values = pal_vert) +
  labs(x = "Theoretical quantiles (standard normal)", y = "Sample quantiles (degree of spondylolisthesis, mm)",
       title = "Degree of spondylolisthesis is far from normal in the abnormal group",
       subtitle = "Line = quantiles of a normal distribution with the same mean and SD; heavy right tail with one extreme patient (418 mm)")
savep(figV2, "vert_fig2_qq_spondylolisthesis", 12, 5.6)

# ---- Fig V3: Spearman correlation heatmap -------------------------------------------
M <- as.matrix(vert[, feats]); R <- cor(M, method = "spearman")
P <- matrix(NA, 6, 6); for (i in 1:6) for (j in 1:6) P[i, j] <- suppressWarnings(cor.test(M[, i], M[, j], method = "spearman", exact = FALSE)$p.value)
cm <- expand.grid(xi = 1:6, yi = 1:6) |> filter(xi < yi) |>
  mutate(rho = R[cbind(xi, yi)], p = P[cbind(xi, yi)],
         x = factor(nice[feats[xi]], levels = nice), y = factor(nice[feats[yi]], levels = rev(nice)),
         label = ifelse(p < 0.05, sprintf("%.2f*", rho), sprintf("%.2f", rho)),
         txt = ifelse(abs(rho) > 0.55, "white", "#0b0b0b"))
write.csv(cm |> select(x, y, rho, p), file.path(res, "vert_spearman.csv"), row.names = FALSE)
figV3 <- ggplot(cm, aes(x, y, fill = rho)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = label, colour = txt), size = 4) + scale_colour_identity() +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", limits = c(-1, 1), name = "Spearman rho") +
  scale_x_discrete(position = "top") + coord_equal() +
  labs(x = NULL, y = NULL, title = "Pelvic incidence is strongly correlated with sacral slope and lumbar lordosis",
       subtitle = "Spearman rank correlations between the six measurements (n = 310); * p < 0.05") +
  theme(axis.text.x.top = element_text(angle = 35, hjust = 0), panel.grid = element_blank())
savep(figV3, "vert_fig3_spearman_heatmap", 9.5, 7)

# ---- Fig V4: scatter pelvic incidence vs sacral slope + regression line ---------------
ct <- cor.test(vert$pelvic_incidence, vert$sacral_slope); rs <- cor.test(vert$pelvic_incidence, vert$sacral_slope, method = "spearman", exact = FALSE)
fit <- lm(sacral_slope ~ pelvic_incidence, data = vert); co <- coef(fit); r2 <- summary(fit)$r.squared
writeLines(c(capture.output(print(ct)), capture.output(print(rs)), capture.output(print(summary(fit)))),
           file.path(res, "vert_regression_pi_ss.txt"))
figV4 <- ggplot(vert, aes(pelvic_incidence, sacral_slope)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.8, fill = "grey80") +
  geom_point(aes(colour = class), size = 1.9, alpha = 0.75) +
  scale_colour_manual(values = pal_vert, name = NULL) +
  labs(x = "Pelvic incidence (degrees)", y = "Sacral slope (degrees)",
       title = "Sacral slope increases linearly with pelvic incidence in both groups",
       subtitle = sprintf("Pearson r = %.2f, Spearman rho = %.2f (both p < 0.001), n = %d\nLeast-squares line: Sacral slope = %.2f + %.2f x Pelvic incidence, R² = %.2f",
                          ct$estimate, rs$estimate, nrow(vert), co[1], co[2], r2),
       caption = "Grey band = 95% confidence interval of the line. Pelvic incidence = pelvic tilt + sacral slope by definition, which explains the strong correlation.")
savep(figV4, "vert_fig4_scatter_pi_ss", 10, 6.5)

# ---- Fig V5: parallel coordinates plot (z-scores) -------------------------------------
z <- vert |> mutate(id = row_number()) |> mutate(across(all_of(feats), ~ as.numeric(scale(.x)))) |>
  pivot_longer(all_of(feats), names_to = "variable") |> mutate(variable = factor(nice[variable], levels = nice))
zmed <- z |> group_by(class, variable) |> summarise(value = median(value), .groups = "drop")
figV5 <- ggplot(z, aes(variable, value, group = id, colour = class)) +
  geom_line(alpha = 0.18, linewidth = 0.4) +
  geom_line(data = zmed, aes(group = class), linewidth = 1.6) +
  geom_point(data = zmed, aes(group = class), size = 2.6) +
  scale_colour_manual(values = pal_vert, name = NULL) +
  labs(x = NULL, y = "Standardised value, z = (x - mean) / SD",
       title = "Biomechanical profiles: abnormal patients sit higher on every measurement except pelvic radius",
       subtitle = "One thin line per patient (n = 310); thick lines = group medians") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
savep(figV5, "vert_fig5_parallel_coordinates", 12, 6)

# ---- Fig V6: histogram of degree of spondylolisthesis ----------------------------------
figV6 <- ggplot(vert, aes(degree_spondylolisthesis, fill = class)) +
  geom_histogram(bins = 45, colour = "white", linewidth = 0.2, position = "identity", alpha = 0.7) +
  scale_fill_manual(values = pal_vert, name = NULL) +
  labs(x = "Degree of spondylolisthesis (mm)", y = "Number of patients",
       title = "Degree of spondylolisthesis: near zero in normal patients, long right tail in abnormal patients",
       subtitle = sprintf("Mean = %.1f, median = %.1f, skewness = %.1f (n = 310)",
                          mean(vert$degree_spondylolisthesis), median(vert$degree_spondylolisthesis), skew(vert$degree_spondylolisthesis)))
savep(figV6, "vert_fig6_histogram_spondylolisthesis", 11, 5.2)

# ---- Liver dataset: Wilcoxon tests for the five tests shown in the deck -----------------
liver <- read.csv(file.path(root, "data/raw/indian_liver_patient.csv")) |>
  mutate(status = factor(Dataset, levels = c(1, 2), labels = c("Liver patient", "No liver disease")))
lv <- c(Total_Bilirubin = "Total bilirubin (mg/dL)", Alkaline_Phosphotase = "Alkaline phosphatase (IU/L)",
        Alamine_Aminotransferase = "ALT (IU/L)", Aspartate_Aminotransferase = "AST (IU/L)", Albumin = "Albumin (g/dL)")
liver_tab <- lapply(names(lv), function(v) { w <- wilcox.test(liver[[v]] ~ liver$status)
  data.frame(test = lv[[v]], median_liver = median(liver[[v]][liver$status == "Liver patient"]),
             median_no = median(liver[[v]][liver$status == "No liver disease"]), W = w$statistic, p = w$p.value) }) |> bind_rows()
write.csv(liver_tab, file.path(res, "liver_group_comparison_wilcoxon.csv"), row.names = FALSE); print(liver_tab)
sk_b <- c(raw = skew(liver$Total_Bilirubin), log10 = skew(log10(liver$Total_Bilirubin)))
write.csv(data.frame(scale = names(sk_b), skewness = sk_b, shapiro_W = c(shapiro.test(liver$Total_Bilirubin)$statistic, shapiro.test(log10(liver$Total_Bilirubin))$statistic)),
          file.path(res, "liver_bilirubin_skewness.csv"), row.names = FALSE)
message("presentation figures done")
