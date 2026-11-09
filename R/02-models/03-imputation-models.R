# Models that estimate missing platelet values within the model (sensitivity analysis;
# Tables S2, S3, S8, S9, S12 and Figure S8)
#
# Platelet predictors enter through brms::mi(): missing values are treated as parameters and
# estimated jointly with the outcome, using Gaussian sub-models with time and a participant random
# intercept. SLOW: the first six models use adapt_delta = 0.999 and max_treedepth = 50 and take
# several hours each; the last four take 20-40 minutes each.
#
#   mi_frail, mi_function                          adjusted models for all visits with the outcome (Tables S12, S8, S9)
#   mi_frail_platelets, mi_function_platelets      + total platelet count (Table S2)
#   mi_frail_logratio, mi_function_logratio        log-ratio of CD41+CD62+ to CD41+CD62- (Table S3)
#   mi_frail_97, mi_function_97                    restricted to participants with >= 1 platelet measurement
#   mi_frail_conf, mi_function_conf                + diabetes + statin use

source(here::here("R/functions.R"))
library(brms)

d <- prepare_analysis_data()

plt_resp <- c("cd41positivecd62negativestd", "cd41positivecd62positivestd")
fn_resp <- c("dinamometeraveragestd", "sftcststd", "sftelbowflexionstd", "sftfutstd")

prior_mi_frail <- function(resp) {
  c(set_prior("normal(0, 1.5)", class = "b", resp = resp),
    set_prior("normal(1, 1.5)", class = "sd", lb = 0, resp = resp),
    set_prior("normal(1, 1.5)", class = "sigma", lb = 0, resp = setdiff(resp, "frailtotal")))
}

submodel <- function(variable) bf(as.formula(paste0(variable, " | mi() ~ time + (1 | id)")), family = gaussian())
neg <- submodel("cd41_positive_cd62_negative_std")
pos <- submodel("cd41_positive_cd62_positive_std")

mi_frail_formula <- function(extra = "", plt_term = "", extra_sub = NULL) {
  f <- bf(as.formula(paste0("frail_total | trials(5) ~ time * (mi(cd41_positive_cd62_negative_std) + ",
                            "mi(cd41_positive_cd62_positive_std) + sex + age_std + imc_std", plt_term, ")", extra,
                            " + (1 | id)")), family = binomial("logit"))
  if (!is.null(extra_sub)) f + extra_sub + neg + pos else f + neg + pos
}
mi_function_formula <- function(extra = "", plt_term = "", extra_sub = NULL) {
  f <- bf(as.formula(paste0("mvbind(dinamometer_average_std, sft_cst_std, sft_elbow_flexion_std, sft_fut_std) ~ ",
                            "time * (mi(cd41_positive_cd62_negative_std) + mi(cd41_positive_cd62_positive_std) + ",
                            "sex + age_std + imc_std", plt_term, ")", extra, " + (1 | id)")), family = gaussian()) +
    set_rescor(TRUE)
  if (!is.null(extra_sub)) f + extra_sub + neg + pos else f + neg + pos
}

slow <- list(iter = 10000, warmup = 5000, adapt_delta = 0.999, max_treedepth = 50)
fast <- list(iter = 6000, warmup = 3000, adapt_delta = 0.99, max_treedepth = 12)
fit_slow <- function(name, formula, data, prior) do.call(fit_model, c(list(name, formula, data, prior), slow))
fit_fast <- function(name, formula, data, prior) do.call(fit_model, c(list(name, formula, data, prior), fast))

# Adjusted models (all visits with the outcome) ---------------------------------------------------------

fit_slow("mi_frail", mi_frail_formula(), d, prior_mi_frail(c("frailtotal", plt_resp)))
fit_slow("mi_function", mi_function_formula(), d, prior_for(c(fn_resp, plt_resp)))

# + total platelet count (Table S2) ---------------------------------------------------------------------------

platelets_sub <- bf(platelets | mi() ~ time + (1 | id), family = gaussian())
fit_slow("mi_frail_platelets",
         mi_frail_formula(plt_term = " + mi(platelets)", extra_sub = platelets_sub), d,
         prior_mi_frail(c("frailtotal", plt_resp, "platelets")))
fit_slow("mi_function_platelets",
         mi_function_formula(plt_term = " + mi(platelets)", extra_sub = platelets_sub), d,
         prior_for(c(fn_resp, plt_resp, "platelets")))

# Log-ratio of CD41+CD62+ to CD41+CD62- (Table S3): 15,000 iterations ------------------------------------------

logratio_sub <- bf(ilr_positive_to_positive_std | mi() ~ time + (1 | id), family = gaussian())
fit_model("mi_frail_logratio",
          bf(frail_total | trials(5) ~ time * (mi(ilr_positive_to_positive_std) + sex + age_std + imc_std) + (1 | id),
             family = binomial("logit")) + logratio_sub,
          d, prior_mi_frail(c("frailtotal", "ilrpositivetopositivestd")),
          iter = 15000, warmup = 5000, adapt_delta = 0.999, max_treedepth = 50)
fit_model("mi_function_logratio",
          bf(mvbind(dinamometer_average_std, sft_cst_std, sft_elbow_flexion_std, sft_fut_std) ~
               time * (mi(ilr_positive_to_positive_std) + sex + age_std + imc_std) + (1 | id), family = gaussian()) +
            set_rescor(TRUE) + logratio_sub,
          d, prior_for(c(fn_resp, "ilrpositivetopositivestd")),
          iter = 15000, warmup = 5000, adapt_delta = 0.999, max_treedepth = 50)

# Participants with at least one platelet measurement -------------------------------------------------------------

ids_platelets <- d[!is.na(cd41_positive_cd62_negative) & !is.na(cd41_positive_cd62_positive), unique(id)]
d97 <- d[id %in% ids_platelets]
fit_fast("mi_frail_97", mi_frail_formula(), d97, prior_mi_frail(c("frailtotal", plt_resp)))
fit_fast("mi_function_97", mi_function_formula(), d97, prior_for(c(fn_resp, plt_resp)))

# + diabetes and statin use ---------------------------------------------------------------------------------------------

fit_fast("mi_frail_conf", mi_frail_formula(extra = " + diabetes + statins"), d, prior_mi_frail(c("frailtotal", plt_resp)))
fit_fast("mi_function_conf", mi_function_formula(extra = " + diabetes + statins"), d, prior_for(c(fn_resp, plt_resp)))
