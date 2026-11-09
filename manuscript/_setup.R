# Shared setup for the manuscript, the supplementary material and the response letter ------------
#
# Everything reported in the documents comes from the objects created by the scripts in R/
# (output/ and manuscript/figures/). Nothing is re-fitted here.

# Make sure strings (e.g. "±", superscripts) are handled as UTF-8
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))

suppressPackageStartupMessages({
  library(data.table)
  library(flextable)
})

proj <- here::here()
source(file.path(proj, "R/functions.R"))

tables_dir <- file.path(proj, "output/tables")
table_file <- function(...) file.path(tables_dir, ...)

# Data ----------------------------------------------------------------------------------------------------

load(file.path(proj, "data/platelets_d.RData"))

# Spearman correlations between age and platelets (Figure S3 caption)
age_cor <- age_platelet_correlations(platelets_d)$correlations
rho_p <- function(subtype) {
  r <- age_cor[age_cor$`Platelet subtype` == subtype]
  sprintf("rho = %.2f, p-value = %.3f", round(r$rho, 2) + 0, r$p.value)
}

# Prevalence of a condition ("Yes") among participants, from the first record of each participant (as in Table 1)
prevalence <- function(variable) {
  first <- platelets_d[, .SD[1], id]
  sprintf("%.1f%%", 100 * mean(first[[variable]] == "Yes"))
}

# Posterior summaries of every model (output/estimates/model-estimates.csv) ---------------------------------

estimates <- read_estimates()
posterior <- split(estimates, by = "model")

# Inline statistics -----------------------------------------------------------------------------------------------
#
# Formats a posterior summary as used in the Results section:
#   "<beta> = est, CI95% [lo, hi], pd = x%, ps = y%"
#
#   var    parameter name in the model-estimate store
#   model  model key (e.g. "main_frail", "main_function", "plt_neg")
#   sym    "beta" (default) or "r" (residual correlations)
#   pct    "1dp" -> 100.0%; "short" -> 100% (only affects values equal to 100)
stat <- function(var, model, sym = "beta", pct = "1dp") {
  post <- posterior[[model]]
  idx <- which(post$var == var)
  stopifnot(length(idx) == 1)
  row <- post[idx]
  ci <- as.numeric(strsplit(gsub("\\[|\\]", "", row$ci), ", ")[[1]])
  p <- function(x) {
    x <- round(x * 100, 1)
    if (pct == "short" && x == 100) "100" else sprintf("%.1f", x)
  }
  symbol <- if (sym == "r") "$r$" else "$\\beta$"
  sprintf("%s = %.2f, CI~95%%~[%s, %s], pd = %s%%, ps = %s%%",
          symbol, row$estimate, as.character(ci[1]), as.character(ci[length(ci)]), p(row$pd), p(row$ps))
}

# "est [lo, hi]" for one parameter
ec <- function(var, model, digits = 2) {
  post <- posterior[[model]]
  idx <- which(post$var == var)
  stopifnot(length(idx) == 1)
  row <- post[idx]
  ci <- as.numeric(strsplit(gsub("\\[|\\]", "", row$ci), ", ")[[1]])
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), row$estimate, ci[1], ci[length(ci)])
}

# Probability of direction (as a percentage) of one parameter
pd_cell <- function(var, model) {
  post <- posterior[[model]]
  idx <- which(post$var == var)
  stopifnot(length(idx) == 1)
  sprintf("%.0f%%", 100 * post$pd[idx])
}

# Tables ----------------------------------------------------------------------------------------------------------------

# Light-grey ruled table (as in the submitted documents)
#   widths: column widths in inches; align: one of "left"/"center"/"right" per column
style_table <- function(ft, widths, align, font = "Calibri", size = 8, color = "#333333",
                        rule = "#D3D3D3", frame = "#A8A8A8", header_rows = 1) {
  ft <- flextable::font(ft, fontname = font, part = "all")
  ft <- flextable::fontsize(ft, size = size, part = "all")
  ft <- flextable::color(ft, color = color, part = "all")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::padding(ft, padding.top = 2, padding.bottom = 2, padding.left = 3, padding.right = 3, part = "all")
  ft <- flextable::border_remove(ft)
  ft <- flextable::hline(ft, border = officer::fp_border(color = rule, width = 0.75), part = "all")
  ft <- flextable::hline_top(ft, border = officer::fp_border(color = frame, width = 1.5), part = "all")
  ft <- flextable::hline_bottom(ft, border = officer::fp_border(color = frame, width = 1.5), part = "body")
  for (j in seq_along(align)) ft <- flextable::align(ft, j = j, align = align[j], part = "all")
  ft <- flextable::width(ft, width = widths)
  ft <- flextable::set_table_properties(ft, layout = "fixed", align = "left")
  flextable::line_spacing(ft, space = 1, part = "all")
}

# Cell of a sensitivity comparison table by outcome, term and model column
#   name: "s8" (frailty) or "s9" (physical function)
sens_cell <- function(name, outcome, term, model) {
  file <- c(s8 = "models/s8-frailty.csv", s9 = "models/s9-function.csv")[[name]]
  tab <- fread(table_file(file), na.strings = "")
  tab[Outcome == outcome & Term == term][[model]]
}
# Effect in original units (Table S9b)
unit_cell <- function(model, outcome, platelet) {
  tab <- fread(table_file("models/s9b-original-units.csv"))
  r <- tab[Model == model & Outcome == outcome & Platelet == platelet]
  sprintf("%s %s", r$`Effect per 1 SD (original units)`, r$Unit)
}
# SD of a platelet percentage (percentage points), from the original-units table
sd_cell <- function(platelet) {
  tab <- fread(table_file("models/s9b-original-units.csv"))
  tab[Platelet == platelet][[5]][1]
}
# Value from the gate-composition table by visit (Table S11)
gate_cell <- function(measure, column) {
  tab <- fread(table_file("gate/s11-by-visit.csv"))
  tab[Measure == measure][[column]]
}
# Assessment months and gate summaries by visit (Table S11b)
date_cell <- function(visit_name, column) {
  tab <- fread(table_file("gate/s11-dates.csv"))
  tab[tab$Visit == visit_name][[column]]
}
# CD41+ share by month of acquisition (Table S11c)
month_cell <- function(visit_name, month) {
  tab <- fread(table_file("gate/s11-months.csv"))
  r <- tab[tab$Visit == visit_name & tab$Month == month]
  sub(" \\[.*", "", r[[4]])
}
# Discrete posterior predictive counts (Figure S6)
ppc_cell <- function(model_name, stratum_name, score) {
  tab <- fread(table_file("diagnostics/ppc-frail-counts.csv"))
  idx <- which(tab$model == model_name & tab$stratum == stratum_name & tab$score_label == score)
  stopifnot(length(idx) == 1)
  r <- tab[idx]
  sprintf("observed %d, predicted %.1f [95%% interval %.0f to %.0f]", r$observed, r$predicted, r$lower, r$upper)
}
