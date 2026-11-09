# Comparison of the platelet coefficients across models (Tables S8, S8b, S9, S9b, S10)
#
#   S8   frailty: main (complete-case) model, alternative models, + confounders, imputation-based models
#   S8b  frailty: odds ratios per 1 SD of each platelet percentage
#   S9   physical function: main model, + confounders, imputation-based models
#   S9b  physical function: effects per 1 SD in original units
#   S10  log-ratio of CD41+CD62+ to CD41+CD62- events (complete cases)

source(here::here("R/functions.R"))

out_dir <- here("output/tables/models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Helpers ---------------------------------------------------------------------------------------------

# Harmonize parameter names across models fitted with / without mi() -> response, term
tidy_est <- function(model) {
  x <- read_estimates(model)
  x[, term_raw := sub("^(bsp|b)_", "", var)]
  x[, term_raw := gsub("mi(?=cd41|ilr)", "", term_raw, perl = TRUE)]
  x[, response := fcase(
    grepl("^frailtotal_|^frail3cat_", term_raw), sub("_.*", "", term_raw),
    grepl("^(dinamometeraveragestd|sftcststd|sftelbowflexionstd|sftfutstd)_", term_raw), sub("_.*", "", term_raw),
    default = NA_character_)]
  # univariate frailty models have no response prefix
  x[is.na(response) & !grepl("^rescor|^cd41positive|^ilrpositive", term_raw) &
      grepl("^(Intercept|time|sex|age|imc|diabetes|statins|cd41_|ilr_)", term_raw), response := "frailtotal"]
  x[, term := sub("^(frailtotal|frail3cat|dinamometeraveragestd|sftcststd|sftelbowflexionstd|sftfutstd)_", "", term_raw)]
  x[, ci := collapse_hdi(ci)]
  x[, list(response, term, estimate, ci, pd)]
}

resp_label <- c(frailtotal = "FRAIL score", frail3cat = "FRAIL score",
                dinamometeraveragestd = "Handgrip strength", sftcststd = "Chair sit-to-stand",
                sftelbowflexionstd = "Arm curl", sftfutstd = "Timed up-and-go")

term_label <- function(term) {
  out <- sub("cd41_positive_cd62_negative_std", "CD41+CD62-", term)
  out <- sub("cd41_positive_cd62_positive_std", "CD41+CD62+", out)
  out <- sub("ilr_positive_to_positive_std", "log-ratio CD41+CD62+/CD41+CD62-", out)
  out <- sub("^time.L:(.*)$", "Time [Linear] x \\1", out)
  out <- sub("^(.*):time.L$", "Time [Linear] x \\1", out)
  out <- sub("^time.Q:(.*)$", "Time [Quadratic] x \\1", out)
  sub("^(.*):time.Q$", "Time [Quadratic] x \\1", out)
}

platelet_pair <- c("cd41_positive_cd62_negative_std", "cd41_positive_cd62_positive_std")

key_terms <- function(x, platelet = platelet_pair) {
  pat <- paste0("^(", paste(platelet, collapse = "|"), ")(:time\\.[LQ])?$|^time\\.[LQ]:(", paste(platelet, collapse = "|"), ")$")
  x[grepl(pat, term)]
}

# models: named character vector (column label = model key in the estimate store)
compare <- function(models, platelet = platelet_pair) {
  long <- rbindlist(lapply(names(models), function(label) {
    t <- key_terms(tidy_est(models[[label]]), platelet)
    t[, model := label][]
  }))
  long[, cell := sprintf("%.2f %s", estimate, ci)]
  long[, `:=`(Outcome = resp_label[response], Term = term_label(term))]
  wide <- dcast(long, Outcome + Term ~ model, value.var = "cell")
  setcolorder(wide, c("Outcome", "Term", names(models)))
  wide[]
}

# Tables S8, S9, S10 -------------------------------------------------------------------------------------

fwrite(compare(c("Main (complete cases)" = "main_frail",
                 "Beta-binomial" = "cc_frail_bb",
                 "Ordinal (0, 1, >=2)" = "cc_frail_ord",
                 "+ diabetes, statins" = "cc_frail_conf",
                 "Imputation-based" = "mi_frail",
                 "Imputation-based, participants with platelet data" = "mi_frail_97")),
       file.path(out_dir, "s8-frailty.csv"))

fwrite(compare(c("Main (complete cases)" = "main_function",
                 "+ diabetes, statins" = "cc_function_conf",
                 "Imputation-based" = "mi_function",
                 "Imputation-based, participants with platelet data" = "mi_function_97",
                 "Imputation-based + diabetes, statins" = "mi_function_conf")),
       file.path(out_dir, "s9-function.csv"))

logratio <- "ilr_positive_to_positive_std"
fwrite(rbind(compare(c("Frailty (log-ratio, complete cases)" = "cc_frail_logratio"), logratio),
             compare(c("Physical function (log-ratio, complete cases)" = "cc_function_logratio"), logratio),
             fill = TRUE),
       file.path(out_dir, "s10-logratio.csv"))

# Effects per 1 SD in original units (Table S9b) and odds ratios (Table S8b) ----------------------------------
# Standardized coefficients x SD of the outcome (median and limits of the 95% HDI are linear in the
# coefficient, so the back-transformation is exact).

load(here("data/platelets_d.RData"))
sd_outcome <- c(dinamometeraveragestd = sd(platelets_d$dinamometer_average, na.rm = TRUE),
                sftcststd = sd(platelets_d$sft_cst, na.rm = TRUE),
                sftelbowflexionstd = sd(platelets_d$sft_elbow_flexion, na.rm = TRUE),
                sftfutstd = sd(platelets_d$sft_fut, na.rm = TRUE))
unit_label <- c(dinamometeraveragestd = "kg", sftcststd = "repetitions", sftelbowflexionstd = "repetitions", sftfutstd = "seconds")
sd_platelet <- c("CD41+CD62-" = sd(platelets_d$cd41_positive_cd62_negative, na.rm = TRUE),
                 "CD41+CD62+" = sd(platelets_d$cd41_positive_cd62_positive, na.rm = TRUE))

limits <- function(t) {
  t[, c("lo", "hi") := list(as.numeric(sub("\\[(.*), .*", "\\1", ci)), as.numeric(sub(".*, (.*)\\]", "\\1", ci)))]
}

to_units <- function(model, label) {
  t <- limits(key_terms(tidy_est(model))[!grepl("time", term) & response %in% names(sd_outcome)])
  t[, `:=`(Model = label, Outcome = resp_label[response], Unit = unit_label[response],
           Platelet = term_label(term), sdy = sd_outcome[response])]
  t[, `:=`(Platelet_SD = sprintf("%.1f", sd_platelet[Platelet]),
           Effect = sprintf("%.2f [%.2f, %.2f]", estimate * sdy, lo * sdy, hi * sdy))]
  t[, list(Model, Outcome, Unit, Platelet, `SD of platelet proportion (%)` = Platelet_SD,
           `Effect per 1 SD (original units)` = Effect)]
}
fwrite(rbind(to_units("main_function", "Main (complete cases)"),
             to_units("cc_function_conf", "+ diabetes, statins"),
             to_units("mi_function", "Imputation-based")),
       file.path(out_dir, "s9b-original-units.csv"))

to_or <- function(model, label) {
  t <- limits(key_terms(tidy_est(model))[!grepl("time", term) & response == "frailtotal"])
  t[, list(Model = label, Platelet = term_label(term),
           `Odds ratio per 1 SD` = sprintf("%.2f [%.2f, %.2f]", exp(estimate), exp(lo), exp(hi)))]
}
fwrite(rbind(to_or("main_frail", "Main (complete cases)"),
             to_or("cc_frail_bb", "Beta-binomial"),
             to_or("cc_frail_conf", "+ diabetes, statins"),
             to_or("mi_frail", "Imputation-based")),
       file.path(out_dir, "s8b-odds-ratios.csv"))
