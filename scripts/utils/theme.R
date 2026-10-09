# ---------------------------------------------------------------------------
# Shared plotting theme, colour palette and figure-saving helpers
# BIOL 520/820 - Assignment 1 (Data Visualization)
#
# Colour system: a colour-blind-checked categorical palette (8 fixed slots),
# a single-hue sequential ramp (blue) and a blue <-> grey <-> red diverging
# scale. Hues are assigned to categories in a FIXED order so that the same
# category always has the same colour across every figure in the report.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# ---- ink / chrome ----------------------------------------------------------
ink <- list(
  primary   = "#0b0b0b",
  secondary = "#52514e",
  muted     = "#898781",
  grid      = "#e1e0d9",
  axis      = "#c3c2b7",
  surface   = "#ffffff"
)

# ---- categorical palette (fixed slot order) --------------------------------
pal_cat <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a",
             yellow = "#eda100", magenta = "#e87ba4", green = "#008300",
             violet = "#4a3aa7", red = "#e34948")

# Dataset 1 - breast cancer subtypes (biological display order).
subtype_levels <- c("normal", "luminal_A", "luminal_B", "HER", "basal", "cell_line")
subtype_labels <- c("Normal", "Luminal A", "Luminal B", "HER2-enriched",
                    "Basal-like", "Cell line")
pal_subtype <- c("Normal"        = pal_cat[["green"]],
                 "Luminal A"     = pal_cat[["blue"]],
                 "Luminal B"     = pal_cat[["orange"]],
                 "HER2-enriched" = pal_cat[["aqua"]],
                 "Basal-like"    = pal_cat[["magenta"]],
                 "Cell line"     = pal_cat[["yellow"]])
# Point shapes used as a second (non-colour) identity channel in scatter plots.
shape_subtype <- c("Normal" = 16, "Luminal A" = 17, "Luminal B" = 15,
                   "HER2-enriched" = 18, "Basal-like" = 8, "Cell line" = 4)

# Dataset 2 - survival outcome (two groups).
pal_outcome <- c("Survived" = pal_cat[["blue"]], "Died" = pal_cat[["orange"]])
shape_outcome <- c("Survived" = 16, "Died" = 17)

# Ordered groups (e.g. ejection-fraction bands): one hue, light -> dark.
pal_ordinal3 <- c("#86b6ef", "#2a78d6", "#0d366b")

# Sequential (magnitude) and diverging (polarity) scales.
pal_seq  <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
pal_div  <- c(low = "#2a78d6", mid = "#f0efec", high = "#e34948")
col_ns   <- "#c3c2b7"   # "not significant" / de-emphasised marks

# ---- ggplot2 theme ---------------------------------------------------------
theme_biol820 <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      text                = element_text(colour = ink$primary),
      plot.background     = element_rect(fill = ink$surface, colour = NA),
      panel.background    = element_rect(fill = ink$surface, colour = NA),
      panel.grid.major    = element_line(colour = ink$grid, linewidth = 0.3),
      panel.grid.minor    = element_blank(),
      axis.ticks          = element_blank(),
      axis.text           = element_text(colour = ink$secondary, size = rel(0.85)),
      axis.title          = element_text(colour = ink$secondary, size = rel(0.9)),
      plot.title          = element_text(face = "bold", size = rel(1.15), colour = ink$primary),
      plot.subtitle       = element_text(colour = ink$secondary, size = rel(0.9),
                                         margin = margin(b = 8)),
      plot.caption        = element_text(colour = ink$muted, size = rel(0.75), hjust = 0,
                                         margin = margin(t = 10)),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      legend.position     = "top",
      legend.title        = element_text(size = rel(0.85), colour = ink$secondary),
      legend.text         = element_text(size = rel(0.85), colour = ink$primary),
      legend.key.size     = unit(0.9, "lines"),
      strip.text          = element_text(face = "bold", colour = ink$primary, size = rel(0.9),
                                         hjust = 0),
      strip.background    = element_blank(),
      plot.margin         = margin(12, 14, 10, 12)
    )
}
theme_set(theme_biol820())

# ---- saving helpers --------------------------------------------------------
# Every figure is written twice: a vector PDF (for the report) and a 300-dpi PNG.
fig_dir <- function(root) list(pdf = file.path(root, "figures", "pdf"),
                               png = file.path(root, "figures", "png"))

save_fig <- function(plot, name, root, width = 9, height = 6, dpi = 300) {
  d <- fig_dir(root)
  dir.create(d$pdf, recursive = TRUE, showWarnings = FALSE)
  dir.create(d$png, recursive = TRUE, showWarnings = FALSE)
  pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(d$pdf, paste0(name, ".pdf")), plot, device = pdf_dev,
         width = width, height = height, units = "in", bg = ink$surface)
  png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else grDevices::png
  ggsave(file.path(d$png, paste0(name, ".png")), plot, device = png_dev,
         width = width, height = height, units = "in", dpi = dpi, bg = ink$surface)
  message("saved ", name, " (pdf + png)")
  invisible(file.path(d$pdf, paste0(name, ".pdf")))
}

# For grid objects (e.g. the gtable returned by pheatmap()).
save_grob <- function(grob, name, root, width = 9, height = 6, dpi = 300) {
  d <- fig_dir(root)
  dir.create(d$pdf, recursive = TRUE, showWarnings = FALSE)
  dir.create(d$png, recursive = TRUE, showWarnings = FALSE)
  draw <- function() { grid::grid.newpage(); grid::grid.draw(grob) }
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(file.path(d$pdf, paste0(name, ".pdf")), width = width, height = height, bg = ink$surface)
  } else {
    grDevices::pdf(file.path(d$pdf, paste0(name, ".pdf")), width = width, height = height, bg = ink$surface)
  }
  draw(); dev.off()
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(file.path(d$png, paste0(name, ".png")), width = width, height = height,
                  units = "in", res = dpi, background = ink$surface)
  } else {
    grDevices::png(file.path(d$png, paste0(name, ".png")), width = width, height = height,
                   units = "in", res = dpi, bg = ink$surface)
  }
  draw(); dev.off()
  message("saved ", name, " (pdf + png)")
  invisible(file.path(d$pdf, paste0(name, ".pdf")))
}

# Format p-values for annotations.
fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, format = "f", digits = 3)))

# Locate the project root when a script is run with Rscript (falls back to getwd()).
project_root <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(getwd())
}
