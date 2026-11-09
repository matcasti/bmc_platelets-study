# Sampling diagnostics of the fitted models (Tables S6 and S10b)
#
# Maximum R-hat, minimum bulk and tail ESS, divergent transitions, iterations hitting the maximum
# tree depth, E-BFMI per chain, observations and participants of each model.

source(here::here("R/functions.R"))
library(brms)

diag_row <- function(name, label) {
  fit <- readRDS(here("models", paste0(name, ".rds")))
  s <- summary(fit)
  parameters <- rbind(s$fixed, s$random[[1]], s$spec_pars, s$rescor_pars)
  sf <- fit$fit
  n_iter <- (sf@sim$iter - sf@sim$warmup) * sf@sim$chains
  max_td <- sf@stan_args[[1]]$control$max_treedepth
  data.table(Model = label,
             `Max R-hat` = sprintf("%.3f", max(parameters$Rhat, na.rm = TRUE)),
             `Min bulk ESS` = sprintf("%.0f", min(parameters$Bulk_ESS, na.rm = TRUE)),
             `Min tail ESS` = sprintf("%.0f", min(parameters$Tail_ESS, na.rm = TRUE)),
             `Divergent transitions` = sprintf("%d of %d", rstan::get_num_divergent(sf), n_iter),
             `Max-treedepth hits` = sprintf("%d of %d", rstan::get_num_max_treedepth(sf), n_iter),
             `E-BFMI (per chain)` = paste(sprintf("%.2f", rstan::get_bfmi(sf)), collapse = ", "),
             N = nrow(fit$data), Participants = uniqueN(fit$data$id))
}

models <- c(
  main_frail = "Frailty, main (complete cases)",
  main_function = "Physical function, main (complete cases)",
  cc_frail_bb = "Frailty, beta-binomial",
  cc_frail_ord = "Frailty, collapsed ordinal",
  cc_frail_conf = "Frailty, + diabetes and statins",
  cc_function_conf = "Physical function, + diabetes and statins",
  cc_frail_logratio = "Frailty, log-ratio",
  cc_function_logratio = "Physical function, log-ratio",
  plt_neg = "CD41+CD62- across visits",
  plt_pos = "CD41+CD62+ across visits",
  mi_frail_97 = "Frailty, imputation-based, participants with platelet data",
  mi_function_97 = "Physical function, imputation-based, participants with platelet data",
  mi_frail_conf = "Frailty, imputation-based + diabetes and statins",
  mi_function_conf = "Physical function, imputation-based + diabetes and statins"
)

missing <- names(models)[!file.exists(here("models", paste0(names(models), ".rds")))]
if (length(missing) > 0) {
  message("Fitted models not found in models/ (", paste(missing, collapse = ", "),
          "): the diagnostics table distributed in output/tables/diagnostics/ is kept.")
} else {
  diagnostics <- rbindlist(lapply(names(models), function(m) diag_row(m, models[[m]])))
  dir.create(here("output/tables/diagnostics"), recursive = TRUE, showWarnings = FALSE)
  fwrite(diagnostics, here("output/tables/diagnostics/s6-s10b-diagnostics.csv"))
  diagnostics
}
