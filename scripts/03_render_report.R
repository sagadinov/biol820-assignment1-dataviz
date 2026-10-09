#!/usr/bin/env Rscript
# Knit report/report.Rmd to PDF (needs LaTeX, e.g. TinyTeX) and to Word (needs only pandoc).
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- if (length(this_file)) normalizePath(file.path(dirname(this_file[1]), "..")) else getwd()
extra <- c(path.expand("~/.local/bin"), Sys.glob(path.expand("~/.TinyTeX/bin/*")))
extra <- extra[dir.exists(extra)]
if (length(extra)) Sys.setenv(PATH = paste(c(extra, Sys.getenv("PATH")), collapse = .Platform$path.sep))
options(bitmapType = "cairo")
rmd <- file.path(root, "report", "report.Rmd")
rmarkdown::render(rmd, output_format = "word_document", output_file = "Group_X.docx",
                  output_dir = file.path(root, "report"), knit_root_dir = root, quiet = TRUE)
rmarkdown::render(rmd, output_format = "pdf_document", output_file = "Group_X.pdf",
                  output_dir = file.path(root, "report"), knit_root_dir = root, quiet = TRUE)
message("Report written to report/Group_X.docx and report/Group_X.pdf")
