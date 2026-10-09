#!/usr/bin/env Rscript
# Render the assignment report (report/report.Rmd) to PDF.
# Requires pandoc and a LaTeX distribution (TinyTeX is enough). The script adds
# ~/.local/bin (pandoc) and ~/.TinyTeX/bin/* to PATH when they exist.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- if (length(this_file)) normalizePath(file.path(dirname(this_file[1]), "..")) else normalizePath(getwd())
extra <- c(path.expand("~/.local/bin"), Sys.glob(path.expand("~/.TinyTeX/bin/*")))
extra <- extra[dir.exists(extra)]
if (length(extra)) Sys.setenv(PATH = paste(c(extra, Sys.getenv("PATH")), collapse = .Platform$path.sep))
options(bitmapType = "cairo")   # no X11 on a cluster login node
stopifnot(rmarkdown::pandoc_available("2.11"))
out <- rmarkdown::render(file.path(root, "report", "report.Rmd"),
                         output_file = "Group_X.pdf", output_dir = file.path(root, "report"),
                         knit_root_dir = root, quiet = TRUE)
message("Report written to ", out)
