# Shared colours, ggplot2 theme and a helper that saves every figure as PDF + PNG.
library(ggplot2)

# One fixed colour per group, used in every figure (colour-blind friendly set).
pal_subtype <- c("Normal" = "#008300", "Luminal A" = "#2a78d6", "Luminal B" = "#eb6834",
                 "HER2-enriched" = "#1baf7a", "Basal-like" = "#e87ba4", "Cell line" = "#eda100")
pal_liver   <- c("Liver patient" = "#eb6834", "No liver disease" = "#2a78d6")
grey_mark   <- "#898781"

theme_set(
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "#52514e"),
          plot.caption = element_text(colour = grey_mark, hjust = 0),
          plot.title.position = "plot",
          strip.text = element_text(face = "bold", hjust = 0),
          legend.position = "top")
)

# Save a ggplot as figures/pdf/<name>.pdf and figures/png/<name>.png (300 dpi).
save_fig <- function(p, name, root, width = 8, height = 5) {
  dir.create(file.path(root, "figures", "pdf"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "figures", "png"), recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(root, "figures", "pdf", paste0(name, ".pdf")), p, width = width, height = height,
         device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf, bg = "white")
  ggsave(file.path(root, "figures", "png", paste0(name, ".png")), p, width = width, height = height,
         dpi = 300, bg = "white",
         device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else grDevices::png)
  message("saved ", name)
}

# Project root = parent folder of the scripts/ folder (works with Rscript and in RStudio).
project_root <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(getwd())
}

# "p < 0.001" / "p = 0.034"
fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
