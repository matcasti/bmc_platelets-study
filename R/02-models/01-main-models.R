# Main models (Table S1, Figures 1 and 2)
#
# Bayesian multilevel models fitted to the visits with observed outcomes, covariates and platelet
# percentages (complete cases):
#   main_frail      FRAIL score (binomial, 5 trials) ~ time x (platelets + sex + age + BMI) + (1 | id)
#   main_function   handgrip, chair sit-to-stand, arm curl and timed up-and-go (multivariate Gaussian)
#   plt_neg, plt_pos  each platelet percentage over visits (Gaussian, observed values)
# Sampler: 4 chains x 6,000 iterations (3,000 warmup), adapt_delta = 0.99, max_treedepth = 12.

source(here::here("R/functions.R"))
library(brms)

d <- prepare_analysis_data()
cc <- complete_case_data(d)

fn_resp <- c("dinamometeraveragestd", "sftcststd", "sftelbowflexionstd", "sftfutstd")
prior_binomial <- c(set_prior("normal(0, 1.5)", class = "b"),
                    set_prior("normal(1, 1.5)", class = "sd", lb = 0),
                    set_prior("normal(0, 2.5)", class = "Intercept"))
prior_gaussian_time <- c(set_prior("normal(0, 1.5)", class = "b"),
                         set_prior("normal(1, 1.5)", class = "sd", lb = 0),
                         set_prior("normal(1, 1.5)", class = "sigma", lb = 0))

fit_model(
  "main_frail",
  bf(frail_total | trials(5) ~ time * (cd41_positive_cd62_negative_std + cd41_positive_cd62_positive_std +
                                         sex + age_std + imc_std) + (1 | id), family = binomial("logit")),
  cc$frail, prior_binomial
)

fit_model(
  "main_function",
  bf(mvbind(dinamometer_average_std, sft_cst_std, sft_elbow_flexion_std, sft_fut_std) ~
       time * (cd41_positive_cd62_negative_std + cd41_positive_cd62_positive_std +
                 sex + age_std + imc_std) + (1 | id), family = gaussian()) + set_rescor(TRUE),
  cc$`function`, prior_for(fn_resp)
)

fit_model(
  "plt_neg",
  bf(cd41_positive_cd62_negative_std ~ time + (1 | id), family = gaussian()),
  d[!is.na(cd41_positive_cd62_negative_std)], prior_gaussian_time
)

fit_model(
  "plt_pos",
  bf(cd41_positive_cd62_positive_std ~ time + (1 | id), family = gaussian()),
  d[!is.na(cd41_positive_cd62_positive_std)], prior_gaussian_time
)
