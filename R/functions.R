# Shared functions for the analysis pipeline -------------------------------------------------
#
# All scripts are run from anywhere inside the project: paths are built with here::here().

library(data.table)
library(here)

# Project paths ---------------------------------------------------------------------------------

ESTIMATES_FILE <- here("output/estimates/model-estimates.csv")

# Analysis data ---------------------------------------------------------------------------------

# Loads data/platelets_d.RData and adds the standardized variables used by the models
prepare_analysis_data <- function() {
  load(here("data/platelets_d.RData"))
  d <- copy(platelets_d)

  d[, cd41_positive_cd62_negative_std := scale(cd41_positive_cd62_negative)]
  d[, cd41_positive_cd62_positive_std := scale(cd41_positive_cd62_positive)]
  d[, ilr_positive_to_positive_std := scale(log(cd41_positive_cd62_positive / cd41_positive_cd62_negative))]
  d[, platelets := scale(plt_10_9_l)]
  d[, age_std := scale(age)]
  d[, imc_std := scale(imc)]
  d[, dinamometer_average_std := scale(dinamometer_average)]
  d[, sft_cst_std := scale(sft_cst)]
  d[, sft_elbow_flexion_std := scale(sft_elbow_flexion)]
  d[, sft_fut_std := scale(sft_fut)]

  # FRAIL collapsed to three ordered categories (0, 1, >=2): only 5 observations scored >= 3
  d[, frail_3cat := ordered(pmin(frail_total, 2), levels = 0:2, labels = c("0", "1", ">=2"))]
  d[]
}

# Visits with the outcome(s), covariates and platelet predictors observed (complete cases)
complete_case_data <- function(d) {
  fn_outcomes <- c("dinamometer_average_std", "sft_cst_std", "sft_elbow_flexion_std", "sft_fut_std")
  covs <- c("sex", "age_std", "imc_std")
  plt <- c("cd41_positive_cd62_negative_std", "cd41_positive_cd62_positive_std")
  cc <- function(vars) d[complete.cases(d[, vars, with = FALSE])]
  list(
    frail = cc(c("frail_total", covs, plt)),
    `function` = cc(c(fn_outcomes, covs, plt)),
    frail_logratio = cc(c("frail_total", covs, "ilr_positive_to_positive_std")),
    function_logratio = cc(c(fn_outcomes, covs, "ilr_positive_to_positive_std"))
  )
}

# Posterior summaries ---------------------------------------------------------------------------

bf_01 <- function(x, prior_sigma) {
  dens <- density(x = as.numeric(x))
  dens_at_0 <- approx(dens$x, dens$y, xout = 0, n = 1000, rule = 2)$y
  dnorm(x = 0, mean = 0, sd = prior_sigma) / dens_at_0
}

# Median, 95% HDI, probability of direction (pd), practical significance (ps, ROPE = +-0.1),
# Savage-Dickey Bayes factor, ESS and R-hat of a vector of posterior draws
summary_posterior <- function(x, prior_sigma = 3) {
  estimate <- round(x = median(x), digits = 2)
  ci <- round(x = tidybayes::hdi(x), digits = 2)
  ci <- paste0("[", paste0(ci, collapse = ", "), "]")
  pd_side <- if (estimate < 0) {x < 0} else {x > 0}
  pd <- round(x = sum(pd_side) / length(x), digits = 3)
  ps_side <- if (estimate < 0) {x < -0.1} else {x > 0.1}
  ps <- round(x = sum(ps_side) / length(x), digits = 3)
  bf <- round(x = bf_01(x, prior_sigma), digits = 2)
  ess <- round(x = posterior::ess_basic(x), digits = 1)
  rhat <- round(x = posterior::rhat(x), digits = 3)

  list(estimate = estimate, ci = ci, pd = pd, ps = ps, bf = bf, ess = ess, rhat = rhat)
}

summary_model <- function(model, variable = NULL) {
  brms::as_draws_df(model, variable = variable, regex = TRUE) |>
    lapply(summary_posterior) |>
    rbindlist(idcol = "var")
}

# Model-estimate store ----------------------------------------------------------------------------
# One long csv (output/estimates/model-estimates.csv) with the posterior summaries of every model:
# columns model, var, estimate, ci, pd, ps, bf, ess, rhat.

save_estimates <- function(fit, name) {
  est <- summary_model(fit, variable = "^b|^rescor")
  est <- est[!var %like% "^\\."]
  est[, model := name]
  setcolorder(est, "model")
  dir.create(dirname(ESTIMATES_FILE), recursive = TRUE, showWarnings = FALSE)
  keep <- if (file.exists(ESTIMATES_FILE)) fread(ESTIMATES_FILE)[model != name] else NULL
  store <- rbind(keep, est)
  setorder(store, model)                                  # stable: row order does not depend on the order of the runs
  fwrite(store, ESTIMATES_FILE)
  invisible(est)
}

read_estimates <- function(models = NULL) {
  est <- fread(ESTIMATES_FILE)
  if (is.null(models)) est else est[model %in% models]
}

# Model fitting -------------------------------------------------------------------------------------

# Fits (or loads from models/<name>.rds when nothing changed) a brms model and stores its estimates.
# Set the environment variable BRMS_FILE_REFIT=never to load an existing models/<name>.rds as it is.
fit_model <- function(name, formula, data, prior, iter = 6000, warmup = 3000,
                      adapt_delta = 0.99, max_treedepth = 12) {
  library(brms)
  dir.create(here("models"), showWarnings = FALSE)
  message("Model: ", name)
  fit <- brm(formula = formula, data = data, prior = prior,
             chains = 4, cores = 4, iter = iter, warmup = warmup, seed = 1234,
             control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
             file = here("models", paste0(name, ".rds")),
             file_refit = Sys.getenv("BRMS_FILE_REFIT", "on_change"))
  save_estimates(fit, name)
  invisible(fit)
}

# Priors: N(0, 1.5) for coefficients, N(1, 1.5) (lower bound 0) for standard deviations and residual
# standard deviations, LKJ(2) for residual correlations. `resp` are brms response names.
prior_for <- function(resp, rescor = TRUE) {
  library(brms)
  p <- c(set_prior("normal(0, 1.5)", class = "b", resp = resp),
         set_prior("normal(1, 1.5)", class = "sd", resp = resp, lb = 0),
         set_prior("normal(1, 1.5)", class = "sigma", resp = resp, lb = 0))
  if (rescor) p <- c(p, set_prior("lkj(2)", class = "rescor"))
  p
}

# Parameter tables (Tables S1-S3, S12) -----------------------------------------------------------------

# Collapse multi-modal HDIs (e.g. "[-1.86, 0.51, -0.75, 1.54]") into a single interval spanning
# the outermost limits.
collapse_hdi <- function(ci) {
  vapply(strsplit(gsub("\\[|\\]", "", ci), ", "), function(z) {
    paste0("[", z[1], ", ", z[length(z)], "]")
  }, character(1))
}

format_parameter_table <- function(params) {
  params[, ci := collapse_hdi(ci)]
  params <- params[, list(Response, parameter, estimate, ci, pd, ps, bf, ess, rhat)][order(Response, parameter)]
  params[, `:=`(estimate = sprintf("%.2f", estimate), pd = sprintf("%.3f", pd), ps = sprintf("%.3f", ps),
                bf = sprintf("%.2f", bf), ess = sprintf("%.1f", ess), rhat = sprintf("%.3f", rhat))]
  names(params) <- c("Response", "Parameter", "Estimate", "95% CI", "pd", "ps", "BF", "ESS", "R-hat")
  params[]
}

response_labels <- c(frailtotal = "FRAIL score", platelets = "Absolute platet count",
                     cd41positivecd62negativestd = "CD41+CD62-", cd41positivecd62positivestd = "CD41+CD62+",
                     ilrpositivetopositivestd = "ILR [CD41+CD62+ to CD41+CD62- ratio]",
                     dinamometeraveragestd = "Grip strength", sftcststd = "Chair sit-to-stand",
                     sftelbowflexionstd = "Elbow flexion test", sftfutstd = "Timed-up and go")

# Models that estimate missing platelet values within the model (mi() terms)
parameter_table_mi <- function(models) {
  params <- read_estimates(models)
  # platelet sub-model parameters appear in every model: keep the first model listed
  params <- params[order(match(model, models))][!var %like% "rescor", .SD[1], var]

  params[, parameter := fcase(
    grepl("_Intercept", var), "Intercept",
    grepl("_time.L$", var), "Time [Linear]",
    grepl("_time.Q$", var), "Time [Quadratic]",
    grepl("_sexF$", var), "Sex [Female]",
    grepl("_age_std$", var), "Age",
    grepl("_imc_std$", var), "BMI",
    grepl("_miplatelets$", var), "Absolute platet count",
    grepl("_micd41_positive_cd62_negative_std$", var), "CD41+CD62-",
    grepl("_micd41_positive_cd62_positive_std$", var), "CD41+CD62+",
    grepl("_miilr_positive_to_positive_std$", var), "ILR [CD41+CD62+ to CD41+CD62- ratio]",
    grepl("_miplatelets:time.L$", var), "Time [Linear] x Absolute platet count",
    grepl("_micd41_positive_cd62_negative_std:time.L$", var), "Time [Linear] x CD41+CD62-",
    grepl("_micd41_positive_cd62_positive_std:time.L$", var), "Time [Linear] x CD41+CD62+",
    grepl("_miilr_positive_to_positive_std:time.L$", var), "Time [Linear] x ILR [CD41+CD62+ to CD41+CD62- ratio]",
    grepl("_miplatelets:time.Q$", var), "Time [Quadratic] x Absolute platet count",
    grepl("_micd41_positive_cd62_negative_std:time.Q$", var), "Time [Quadratic] x CD41+CD62-",
    grepl("_micd41_positive_cd62_positive_std:time.Q$", var), "Time [Quadratic] x CD41+CD62+",
    grepl("_miilr_positive_to_positive_std:time.Q$", var), "Time [Quadratic] x ILR [CD41+CD62+ to CD41+CD62- ratio]",
    grepl("_time.L:sexF$", var), "Time [Linear] x Sex [Female]",
    grepl("_time.Q:sexF$", var), "Time [Quadratic] x Sex [Female]",
    grepl("_time.L:age_std$", var), "Time [Linear] x Age",
    grepl("_time.Q:age_std$", var), "Time [Quadratic] x Age",
    grepl("_time.L:imc_std$", var), "Time [Linear] x BMI",
    grepl("_time.Q:imc_std$", var), "Time [Quadratic] x BMI"
  )]
  response_key <- sub("^(bsp|b)_([a-z0-9]+)_.*$", "\\2", params$var)
  params[, Response := response_labels[response_key]]
  format_parameter_table(params)
}

# Models fitted to observed values (no mi() terms): parameter names such as
# "b_sftcststd_time.L:cd41_positive_cd62_negative_std"
parameter_table_observed <- function(models) {
  params <- read_estimates(models)
  params <- params[!var %like% "rescor"]

  params[, response := fcase(
    grepl("^b_(dinamometeraveragestd|sftcststd|sftelbowflexionstd|sftfutstd)_", var), sub("^b_([a-z]+)_.*$", "\\1", var),
    default = "frailtotal")]                      # univariate frailty model: no response prefix
  params[, term := sub("^b_(dinamometeraveragestd|sftcststd|sftelbowflexionstd|sftfutstd)_", "", var)]
  params[, term := sub("^b_", "", term)]

  term_label <- function(x) {
    x <- sub("cd41_positive_cd62_negative_std", "CD41+CD62-", x)
    x <- sub("cd41_positive_cd62_positive_std", "CD41+CD62+", x)
    x <- sub("sexF", "Sex [Female]", x)
    x <- sub("age_std", "Age", x)
    x <- sub("imc_std", "BMI", x)
    x <- sub("^time\\.L:(.*)$", "Time [Linear] x \\1", x)
    x <- sub("^time\\.Q:(.*)$", "Time [Quadratic] x \\1", x)
    x <- sub("^time\\.L$", "Time [Linear]", x)
    sub("^time\\.Q$", "Time [Quadratic]", x)
  }
  params[, parameter := term_label(term)]
  params[, Response := response_labels[response]]
  format_parameter_table(params)
}

# Discrete posterior predictive checks -------------------------------------------------------------------

# Observed vs posterior predictive counts of every score level, overall and by stratum.
#   yrep:   draws x observations matrix of replicated scores
#   y:      observed scores
#   levels: score levels to tabulate
#   strata: vector (same length as y) defining the strata (e.g. visit); NULL = overall only
ppc_counts <- function(yrep, y, levels, strata = NULL, probs = c(0.025, 0.975)) {
  one <- function(idx, label) {
    obs <- vapply(levels, function(l) sum(y[idx] == l), numeric(1))
    rep <- vapply(levels, function(l) rowSums(yrep[, idx, drop = FALSE] == l), numeric(nrow(yrep)))
    data.table(stratum = label, score = levels, n_obs = length(idx), observed = obs,
               predicted = colMeans(rep),
               lower = apply(rep, 2, quantile, probs[1]), upper = apply(rep, 2, quantile, probs[2]))
  }
  out <- list(one(seq_along(y), "Overall"))
  if (!is.null(strata)) {
    for (s in sort(unique(strata))) out[[length(out) + 1]] <- one(which(strata == s), as.character(s))
  }
  rbindlist(out)
}

# Spearman correlations between age and platelet measures (Figure S3 caption) --------------------------------

age_platelet_correlations <- function(data) {
  plot_data <- data[, list(id, Age = age,
                           Platelets = plt_10_9_l,
                           `CD41- CD62+` = cd41_negative_cd62_positive / 100,
                           `CD41+ CD62+` = cd41_positive_cd62_positive / 100,
                           `CD41- CD62-` = cd41_negative_cd62_negative / 100,
                           `CD41+ CD62-` = cd41_positive_cd62_negative / 100)] |>
    na.omit()

  plot_data <- melt.data.table(data = plot_data,
                               id.vars = c("id", "Age", "Platelets"),
                               measure.vars = c("CD41- CD62+", "CD41+ CD62+", "CD41- CD62-", "CD41+ CD62-"),
                               value.name = "Proportion",
                               variable.name = "Platelet subtype")

  total <- cor.test(plot_data$Platelets, plot_data$Age, method = "spearman", exact = FALSE)
  subtype <- plot_data[, {
    cor <- cor.test(Proportion, Age, method = "spearman", exact = FALSE)
    list(rho = unname(cor$estimate), p.value = cor$p.value)
  }, `Platelet subtype`]

  list(data = plot_data,
       correlations = rbind(data.table(`Platelet subtype` = "Total", rho = unname(total$estimate), p.value = total$p.value),
                            subtype, use.names = FALSE))
}
