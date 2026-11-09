# Runs the complete analysis from the analysis dataset to the rendered documents ----------------------
#
#   Rscript run-all.R                    # everything (the imputation-based models take many hours)
#   SKIP_IMPUTATION=1 Rscript run-all.R  # do not fit the imputation-based models (hours): their estimates are read
#                                        # from output/estimates/ (distributed with the repository)
#   BRMS_FILE_REFIT=never Rscript run-all.R  # load models already present in models/ without checking for changes
#   RENDER=0 Rscript run-all.R           # analyses only, without rendering the documents
#
# Fitted models are cached in models/ (brms `file =`): a model is only refitted if its data, formula,
# priors or settings change. Estimates are written to output/estimates/ and tables to output/tables/.
# Run from the project root (the folder that contains this file).

options(warn = 1)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

# Compiling Stan code on macOS fails if Homebrew's LLVM comes first in PATH: use Apple's clang
path <- strsplit(Sys.getenv("PATH"), .Platform$path.sep)[[1]]
Sys.setenv(PATH = paste(path[!grepl("/opt/homebrew/opt/llvm", path)], collapse = .Platform$path.sep))

root <- here::here()
setwd(root)

run <- function(script) {
  message("\n==== ", script, " ====")
  t0 <- Sys.time()
  source(file.path(root, script), local = new.env(), echo = FALSE)
  message(sprintf("     done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# 1. Data ----------------------------------------------------------------------------------------------
# (R/01-data/01-prepare-data.R needs the raw database, which is not distributed: data/platelets_d.RData
#  is provided and is the starting point of every other script)
if (file.exists("data/raw/data-raw.csv") && identical(Sys.getenv("PREPARE_DATA"), "1")) {
  run("R/01-data/01-prepare-data.R")
}
for (s in c("02-table-1", "03-observations-and-missing-data", "04-gate-composition", "05-platelet-change-bootstrap")) {
  run(file.path("R/01-data", paste0(s, ".R")))
}

# 2. Models ----------------------------------------------------------------------------------------------
for (s in c("01-main-models", "02-sensitivity-observed-data", "04-conventional-mixed-models")) {
  run(file.path("R/02-models", paste0(s, ".R")))
}
if (!identical(Sys.getenv("SKIP_IMPUTATION"), "1")) run("R/02-models/03-imputation-models.R")   # slow (hours)

# 3. Summaries -----------------------------------------------------------------------------------------------
run("R/03-summaries/01-model-diagnostics.R")
run("R/03-summaries/02-sensitivity-tables.R")

# 4. Figures -----------------------------------------------------------------------------------------------------
for (s in c("main/fig-1-predictions", "main/fig-2-platelets-over-visits",
            "supplementary/fig-S3-S4-descriptive", "supplementary/fig-S5-S6-S8-posterior-checks",
            "supplementary/fig-S7-S9-gate-composition")) {
  run(file.path("R/04-figures", paste0(s, ".R")))
}

# 5. Documents (Quarto) --------------------------------------------------------------------------------------------
if (!identical(Sys.getenv("RENDER"), "0")) {
  message("\n==== Rendering the documents (Quarto) ====")
  status <- system2("quarto", "render")
  if (status != 0) stop("quarto render failed")
}
message("\nAll done.")
