# Install the CRAN packages used by this project into the user library.
pkgs <- c("ggplot2", "dplyr", "tidyr", "data.table", "scales", "ggrepel", "patchwork",
          "pheatmap", "matrixStats", "survival", "RColorBrewer", "rmarkdown", "knitr", "ragg")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
# A LaTeX distribution is needed for the PDF report; if you have none run:
#   tinytex::install_tinytex()
