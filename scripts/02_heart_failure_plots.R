#!/usr/bin/env Rscript
# ===========================================================================
# Dataset 2 - Heart failure clinical records (299 patients, 13 variables)
# BIOL 520/820 - Assignment 1 (Data Visualization)
#
# Kaggle: https://www.kaggle.com/datasets/andrewmvd/heart-failure-clinical-data
# Original study: Ahmad et al. (2017) PLoS ONE; analysed by Chicco & Jurman (2020).
#
# Figures produced (figures/pdf + figures/png):
#   fig2_1_clinical_violins  continuous clinical variables by outcome   (box/violin)
#   fig2_2_scatter_ef_creat  ejection fraction vs serum creatinine       (scatter plot)
#   fig2_3_km_ejection       Kaplan-Meier curves by ejection fraction    (survival curve)
#   fig2_4_corr_heatmap      Spearman correlation matrix of all variables (heatmap)
# Tables produced (results/): summary by outcome, Wilcoxon tests, quadrant mortality,
#   log-rank test, correlation matrix.
# ===========================================================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(survival)
})
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- if (length(this_file)) normalizePath(file.path(dirname(this_file[1]), "..")) else normalizePath(getwd())
source(file.path(root, "scripts", "utils", "theme.R"))
res_dir <- file.path(root, "results"); dir.create(res_dir, showWarnings = FALSE)
set.seed(820)

# ---------------------------------------------------------------------------
# 1. Load and tidy
# ---------------------------------------------------------------------------
hf <- fread(file.path(root, "data/raw/heart_failure_clinical_records_dataset.csv"))
stopifnot(nrow(hf) == 299, ncol(hf) == 13)
hf <- hf |>
  mutate(outcome   = factor(DEATH_EVENT, levels = c(0, 1), labels = c("Survived", "Died")),
         sex_f     = factor(sex, levels = c(0, 1), labels = c("Female", "Male")),
         ef_group  = cut(ejection_fraction, breaks = c(-Inf, 30, 40, Inf),
                         labels = c("EF ≤ 30%", "EF 31–40%", "EF > 40%")))
var_labels <- c(age = "Age (years)",
                creatinine_phosphokinase = "Creatinine phosphokinase\n(mcg/L, log10)",
                ejection_fraction = "Ejection fraction (%)",
                platelets = "Platelets (10³/µL)",
                serum_creatinine = "Serum creatinine\n(mg/dL, log10)",
                serum_sodium = "Serum sodium (mEq/L)",
                time = "Follow-up time (days)")

# ---- descriptive table by outcome ------------------------------------------
cont_vars <- c("age", "creatinine_phosphokinase", "ejection_fraction", "platelets",
               "serum_creatinine", "serum_sodium", "time")
bin_vars  <- c("anaemia", "diabetes", "high_blood_pressure", "sex", "smoking")
desc_cont <- hf |> select(outcome, all_of(cont_vars)) |>
  mutate(platelets = platelets / 1000) |>           # reported as 10^3 per microlitre
  pivot_longer(-outcome, names_to = "variable") |>
  group_by(variable, outcome) |>
  summarise(stat = sprintf("%.1f [%.1f–%.1f]", median(value), quantile(value, .25), quantile(value, .75)),
            .groups = "drop") |>
  pivot_wider(names_from = outcome, values_from = stat) |>
  mutate(type = "median [IQR]")
desc_bin <- hf |> select(outcome, all_of(bin_vars)) |>
  pivot_longer(-outcome, names_to = "variable") |>
  group_by(variable, outcome) |>
  summarise(stat = sprintf("%d (%.0f%%)", sum(value), 100 * mean(value)), .groups = "drop") |>
  pivot_wider(names_from = outcome, values_from = stat) |>
  mutate(type = "n (%)")
desc <- bind_rows(desc_cont, desc_bin)
fwrite(desc, file.path(res_dir, "ds2_summary_by_outcome.csv"))
fwrite(data.frame(outcome = levels(hf$outcome), n = as.integer(table(hf$outcome))),
       file.path(res_dir, "ds2_outcome_counts.csv"))

# ---------------------------------------------------------------------------
# 2. Figure 2.1 - Continuous clinical variables by outcome (violin + box + points)
# ---------------------------------------------------------------------------
long <- hf |> select(outcome, all_of(cont_vars)) |>
  mutate(creatinine_phosphokinase = log10(creatinine_phosphokinase),
         serum_creatinine = log10(serum_creatinine),
         platelets = platelets / 1000) |>
  pivot_longer(-outcome, names_to = "variable") |>
  mutate(variable = factor(var_labels[variable], levels = var_labels))
wt <- long |> group_by(variable) |>
  summarise(p = wilcox.test(value ~ outcome)$p.value, .groups = "drop") |>
  mutate(label = paste("Mann–Whitney", fmt_p(p)))
fwrite(wt, file.path(res_dir, "ds2_wilcoxon_by_outcome.csv"))
n_out <- table(hf$outcome)

fig2_1 <- ggplot(long, aes(outcome, value)) +
  geom_violin(aes(fill = outcome), colour = NA, alpha = 0.35, width = 0.85, trim = TRUE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = ink$surface, colour = ink$secondary,
               linewidth = 0.35) +
  geom_jitter(aes(colour = outcome), width = 0.13, size = 0.8, alpha = 0.55) +
  geom_text(data = wt, aes(x = 0.5, y = Inf, label = label), hjust = 0, vjust = 1.6,
            size = 2.6, colour = ink$secondary, inherit.aes = FALSE) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = pal_outcome, name = NULL) +
  scale_colour_manual(values = pal_outcome, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)), fill = "none") +
  labs(x = NULL, y = NULL,
       title = "Patients who died during follow-up had lower ejection fraction, higher serum creatinine and shorter follow-up",
       subtitle = sprintf("Heart failure clinical records, n = %d (survived %d, died %d). Violin = distribution, box = median and IQR, points = patients",
                          nrow(hf), n_out[["Survived"]], n_out[["Died"]]),
       caption = "Creatinine phosphokinase and serum creatinine are shown on a log10 scale because of their right-skewed distributions. Two-sided Mann–Whitney U tests.") +
  theme(panel.grid.major.x = element_blank(), legend.position = "top")
save_fig(fig2_1, "fig2_1_clinical_violins", root, width = 12, height = 6.8)

# ---------------------------------------------------------------------------
# 3. Figure 2.2 - Ejection fraction vs serum creatinine, with marginal densities
# ---------------------------------------------------------------------------
ef_cut <- 40; cr_cut <- 1.5
quad <- hf |>
  mutate(quadrant = case_when(ejection_fraction <= ef_cut & serum_creatinine > cr_cut ~ "Low EF, high creatinine",
                              ejection_fraction <= ef_cut & serum_creatinine <= cr_cut ~ "Low EF, normal creatinine",
                              ejection_fraction >  ef_cut & serum_creatinine > cr_cut ~ "Preserved EF, high creatinine",
                              TRUE ~ "Preserved EF, normal creatinine")) |>
  group_by(quadrant) |>
  summarise(n = n(), died = sum(DEATH_EVENT), mortality = died / n, .groups = "drop")
fwrite(quad, file.path(res_dir, "ds2_quadrant_mortality.csv"))
quad_lab <- quad |>
  mutate(x = ifelse(grepl("^Low", quadrant), 14, 81),
         y = ifelse(grepl("high creatinine", quadrant), 9.6, 0.52),
         hjust = ifelse(grepl("^Low", quadrant), 0, 1),
         label = sprintf("%s\nn = %d, died %d (%.0f%%)", quadrant, n, died, 100 * mortality))

p_main <- ggplot(hf, aes(ejection_fraction, serum_creatinine)) +
  annotate("rect", xmin = -Inf, xmax = ef_cut, ymin = cr_cut, ymax = Inf,
           fill = pal_cat[["orange"]], alpha = 0.06) +
  geom_vline(xintercept = ef_cut, linetype = "22", colour = ink$muted, linewidth = 0.35) +
  geom_hline(yintercept = cr_cut, linetype = "22", colour = ink$muted, linewidth = 0.35) +
  geom_point(colour = ink$surface, size = 3.2, position = position_jitter(width = 0.6, height = 0, seed = 3)) +
  geom_point(aes(colour = outcome, shape = outcome), size = 2.2, alpha = 0.85,
             position = position_jitter(width = 0.6, height = 0, seed = 3)) +
  geom_text(data = quad_lab, aes(x, y, label = label, hjust = hjust), vjust = 1, size = 2.7,
            colour = ink$primary, lineheight = 0.95) +
  scale_y_log10(breaks = c(0.5, 1, 1.5, 2, 3, 5, 9)) +
  scale_x_continuous(breaks = seq(10, 80, 10), limits = c(12, 82)) +
  scale_colour_manual(values = pal_outcome, name = NULL) +
  scale_shape_manual(values = shape_outcome, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(x = "Ejection fraction (%)", y = "Serum creatinine (mg/dL, log scale)") +
  theme(legend.position = "bottom", legend.direction = "horizontal")
p_top <- ggplot(hf, aes(ejection_fraction, fill = outcome, colour = outcome)) +
  geom_density(alpha = 0.25, linewidth = 0.5, adjust = 1.1) +
  scale_x_continuous(limits = c(12, 82)) +
  scale_fill_manual(values = pal_outcome, guide = "none") +
  scale_colour_manual(values = pal_outcome, guide = "none") +
  labs(x = NULL, y = "Density", title = "Deaths concentrate where ejection fraction is low and serum creatinine is high",
       subtitle = sprintf("n = %d patients; dashed lines mark EF = %d%% and creatinine = %.1f mg/dL; shaded quadrant = both risk factors present", nrow(hf), ef_cut, cr_cut)) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(),
        axis.title.y = element_text(size = 8))
p_right <- ggplot(hf, aes(serum_creatinine, fill = outcome, colour = outcome)) +
  geom_density(alpha = 0.25, linewidth = 0.5, adjust = 1.1) +
  scale_x_log10() + coord_flip() +
  scale_fill_manual(values = pal_outcome, guide = "none") +
  scale_colour_manual(values = pal_outcome, guide = "none") +
  labs(x = NULL, y = "Density") +
  theme(axis.text = element_blank(), panel.grid = element_blank(), axis.title.x = element_text(size = 8))
fig2_2 <- p_top + plot_spacer() + p_main + p_right +
  plot_layout(ncol = 2, widths = c(5, 1), heights = c(1, 5)) +
  plot_annotation(caption = "Points are jittered horizontally by ±0.6% to reduce over-plotting of repeated ejection-fraction values. Marginal panels: kernel density estimates by outcome.",
                  theme = theme(plot.caption = element_text(colour = ink$muted, size = 8, hjust = 0)))
save_fig(fig2_2, "fig2_2_scatter_ef_creat", root, width = 10, height = 8)

# ---------------------------------------------------------------------------
# 4. Figure 2.3 - Kaplan-Meier survival curves by ejection-fraction band
# ---------------------------------------------------------------------------
sf <- survfit(Surv(time, DEATH_EVENT) ~ ef_group, data = hf)
lr <- survdiff(Surv(time, DEATH_EVENT) ~ ef_group, data = hf)
lr_p <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
strata_names <- sub("^ef_group=", "", names(sf$strata))
km <- data.frame(time = sf$time, surv = sf$surv, lower = sf$lower, upper = sf$upper,
                 n.censor = sf$n.censor, n.event = sf$n.event,
                 group = factor(rep(strata_names, sf$strata), levels = levels(hf$ef_group)))
km0 <- data.frame(time = 0, surv = 1, lower = 1, upper = 1, n.censor = 0, n.event = 0,
                  group = factor(levels(hf$ef_group), levels = levels(hf$ef_group)))
km <- bind_rows(km0, km) |> arrange(group, time)
# step-wise coordinates for confidence ribbons
stepify <- function(d) {
  d <- d[order(d$time), ]
  if (nrow(d) < 2) return(d)
  d2 <- d[-1, ]; d2$surv <- d$surv[-nrow(d)]; d2$lower <- d$lower[-nrow(d)]; d2$upper <- d$upper[-nrow(d)]
  out <- rbind(d, d2); out[order(out$time, -out$surv), ]
}
km_step <- km |> group_by(group) |> group_modify(~ stepify(.x)) |> ungroup()
risk_times <- seq(0, 280, by = 40)
rt <- summary(sf, times = risk_times, extend = TRUE)
risk <- data.frame(time = rt$time, n.risk = rt$n.risk,
                   group = factor(sub("^ef_group=", "", rt$strata), levels = levels(hf$ef_group)))
ends <- km |> group_by(group) |> slice_max(time, n = 1, with_ties = FALSE) |> ungroup() |>
  left_join(as.data.frame(table(group = hf$ef_group)), by = "group") |>
  mutate(label = sprintf("%s (n = %d)", group, Freq))
surv_summary <- summary(sf, times = c(90, 180, 270), extend = TRUE)
fwrite(data.frame(group = sub("^ef_group=", "", surv_summary$strata), time = surv_summary$time,
                  surv = surv_summary$surv, lower = surv_summary$lower, upper = surv_summary$upper,
                  n.risk = surv_summary$n.risk),
       file.path(res_dir, "ds2_km_survival_estimates.csv"))
writeLines(c(sprintf("Log-rank test: chi-square = %.2f, df = %d, p = %.3g", lr$chisq, length(lr$n) - 1, lr_p),
             capture.output(print(sf))), file.path(res_dir, "ds2_km_logrank.txt"))
pal_ef <- setNames(pal_ordinal3, levels(hf$ef_group))
lty_ef <- setNames(c("solid", "42", "12"), levels(hf$ef_group))

p_km <- ggplot(km_step, aes(time, surv, colour = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.10, colour = NA) +
  geom_step(data = km, aes(linetype = group), linewidth = 0.9, direction = "hv") +
  geom_point(data = subset(km, n.censor > 0), shape = 3, size = 1.6, stroke = 0.6, show.legend = FALSE) +
  geom_text(data = ends, aes(x = time + 4, y = surv, label = label), hjust = 0, size = 2.9,
            colour = ink$primary, show.legend = FALSE) +
  annotate("text", x = 2, y = 0.08, hjust = 0, size = 3.2, colour = ink$primary,
           label = sprintf("Log-rank test: χ² = %.1f, df = 2, %s", lr$chisq, fmt_p(lr_p))) +
  scale_colour_manual(values = pal_ef, name = "Ejection fraction") +
  scale_fill_manual(values = pal_ef, name = "Ejection fraction") +
  scale_linetype_manual(values = lty_ef, name = "Ejection fraction") +
  scale_x_continuous(breaks = risk_times, limits = c(0, 345), expand = expansion(mult = c(0.02, 0))) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "Survival probability",
       title = "Survival falls steeply with reduced ejection fraction",
       subtitle = "Kaplan–Meier estimates with 95% confidence bands; '+' marks censored patients") +
  theme(legend.position = "top", legend.key.width = unit(2, "lines"))
p_risk <- ggplot(risk, aes(time, group, label = n.risk)) +
  geom_text(size = 2.9, colour = ink$primary) +
  scale_x_continuous(breaks = risk_times, limits = c(0, 345), expand = expansion(mult = c(0.02, 0))) +
  scale_y_discrete(limits = rev(levels(hf$ef_group))) +
  coord_cartesian(clip = "off") +
  labs(x = "Follow-up time (days)", y = NULL, title = "Number at risk") +
  theme(panel.grid = element_blank(), plot.title = element_text(size = rel(0.85), face = "plain",
                                                                 colour = ink$secondary),
        axis.text.y = element_text(colour = ink$primary, size = 8))
fig2_3 <- p_km / p_risk + plot_layout(heights = c(4.2, 1)) +
  plot_annotation(caption = "Ejection fraction bands follow common clinical cut-offs (≤ 30% severely reduced; 31–40% reduced; > 40% mildly reduced or preserved).\nData: Ahmad et al. (2017); Chicco & Jurman (2020).",
                  theme = theme(plot.caption = element_text(colour = ink$muted, size = 8, hjust = 0)))
save_fig(fig2_3, "fig2_3_km_ejection", root, width = 10, height = 7.2)

# ---------------------------------------------------------------------------
# 5. Figure 2.4 - Spearman correlation heatmap of all variables
# ---------------------------------------------------------------------------
nice <- c(age = "Age", anaemia = "Anaemia", creatinine_phosphokinase = "Creatinine phosphokinase",
          diabetes = "Diabetes", ejection_fraction = "Ejection fraction",
          high_blood_pressure = "High blood pressure", platelets = "Platelets",
          serum_creatinine = "Serum creatinine", serum_sodium = "Serum sodium", sex = "Sex (male)",
          smoking = "Smoking", time = "Follow-up time", DEATH_EVENT = "Death event")
ord <- c("age", "sex", "smoking", "anaemia", "diabetes", "high_blood_pressure",
         "creatinine_phosphokinase", "platelets", "serum_sodium", "serum_creatinine",
         "ejection_fraction", "time", "DEATH_EVENT")
M <- as.matrix(hf[, ..ord])
R <- cor(M, method = "spearman")
P <- matrix(NA, ncol(M), ncol(M), dimnames = dimnames(R))
for (i in seq_len(ncol(M))) for (j in seq_len(ncol(M)))
  P[i, j] <- suppressWarnings(cor.test(M[, i], M[, j], method = "spearman", exact = FALSE)$p.value)
cm <- as.data.frame(as.table(R)) |> setNames(c("x", "y", "rho")) |>
  mutate(p = as.vector(P), xi = as.integer(factor(x, levels = ord)), yi = as.integer(factor(y, levels = ord))) |>
  filter(xi < yi) |>
  mutate(x = factor(nice[as.character(x)], levels = nice[ord]),
         y = factor(nice[as.character(y)], levels = rev(nice[ord])),
         label = ifelse(p < 0.05, sprintf("%.2f*", rho), sprintf("%.2f", rho)),
         txt = ifelse(abs(rho) > 0.45, ink$surface, ink$primary))
fwrite(cm |> select(x, y, rho, p) |> arrange(p), file.path(res_dir, "ds2_spearman_correlations.csv"))

fig2_4 <- ggplot(cm, aes(x, y, fill = rho)) +
  geom_tile(colour = ink$surface, linewidth = 1.2) +
  geom_text(aes(label = label, colour = txt), size = 2.7) +
  scale_colour_identity() +
  scale_fill_gradient2(low = pal_div[["low"]], mid = pal_div[["mid"]], high = pal_div[["high"]],
                       limits = c(-1, 1), name = "Spearman ρ") +
  scale_x_discrete(position = "top") +
  coord_equal() +
  labs(x = NULL, y = NULL,
       title = "Clinical variables are only weakly inter-correlated;\ndeath tracks follow-up time, serum creatinine and ejection fraction",
       subtitle = "Spearman rank correlation between all 13 recorded variables (binary variables coded 0/1)",
       caption = "* p < 0.05 (two-sided test of ρ = 0, no multiplicity correction).\nFollow-up time is shorter for patients who died, hence its strong negative association with the death event.") +
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0), panel.grid = element_blank(),
        legend.position = "right", legend.key.height = unit(1.6, "lines"))
save_fig(fig2_4, "fig2_4_corr_heatmap", root, width = 10, height = 8.5)

writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_ds2.txt"))
message("Dataset 2 done.")
