report_dir <- "reports/final"
report_source <- file.path(report_dir, "final_report.Rmd")
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(report_dir)
Sys.setenv(PROJECT_ROOT = normalizePath("../.."))

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to compile the final report.")
}

rmarkdown::render(
  input = "final_report.Rmd",
  output_format = "pdf_document",
  output_file = "final_report.pdf",
  output_dir = ".",
  clean = TRUE,
  envir = new.env(parent = globalenv())
)
