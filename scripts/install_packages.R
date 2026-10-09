# Install the CRAN packages used by this project (run once).
pkgs <- c("ggplot2", "dplyr", "tidyr", "data.table", "rmarkdown", "knitr", "ragg")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
# For the PDF version of the report a LaTeX distribution is needed: tinytex::install_tinytex()
