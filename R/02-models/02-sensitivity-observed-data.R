# Sensitivity analyses fitted to observed data (complete cases) (Tables S8-S10, Figure S6)
#
#   cc_frail_bb / cc_frail_ord   beta-binomial and collapsed (0, 1, >=2) cumulative-logit FRAIL models
#   cc_frail_conf / cc_function_conf   main models + diabetes + statin use (main effects)
#   cc_frail_logratio / cc_function_logratio   log-ratio of CD41+CD62+ to CD41+CD62- events
# Same sampler settings as the main models.

source(here::here("R/functions.R"))
library(brms)

d <- prepare_analysis_data()
cc <- complete_case_data(d)

fn_resp <- c("dinamometeraveragestd", "sftcststd", "sftelbowflexionstd", "sftfutstd")
prior_binomial <- c(set_prior("normal(0, 1.5)", class = "b"),
                    set_prior("normal(1, 1.5)", class = "sd", lb = 0),
                    set_prior("normal(0, 2.5)", class = "Intercept"))
prior_ordinal <- c(set_prior("normal(0, 1.5)", class = "b"), set_prior("normal(1, 1.5)", class = "sd", lb = 0))

platelet_terms <- "cd41_positive_cd62_negative_std + cd41_positive_cd62_positive_std + sex + age_std + imc_std"
fn_outcomes <- "mvbind(dinamometer_average_std, sft_cst_std, sft_elbow_flexion_std, sft_fut_std)"

fit_model(
  "cc_frail_bb",
  bf(as.formula(paste0("frail_total | trials(5) ~ time * (", platelet_terms, ") + (1 | id)")), family = beta_binomial()),
  cc$frail, prior_binomial
)

fit_model(
  "cc_frail_ord",
  bf(as.formula(paste0("frail_3cat ~ time * (", platelet_terms, ") + (1 | id)")), family = cumulative("logit")),
  cc$frail, prior_ordinal
)

fit_model(
  "cc_frail_conf",
  bf(as.formula(paste0("frail_total | trials(5) ~ time * (", platelet_terms, ") + diabetes + statins + (1 | id)")),
     family = binomial("logit")),
  cc$frail, prior_binomial
)

fit_model(
  "cc_function_conf",
  bf(as.formula(paste0(fn_outcomes, " ~ time * (", platelet_terms, ") + diabetes + statins + (1 | id)")),
     family = gaussian()) + set_rescor(TRUE),
  cc$`function`, prior_for(fn_resp)
)

fit_model(
  "cc_frail_logratio",
  bf(frail_total | trials(5) ~ time * (ilr_positive_to_positive_std + sex + age_std + imc_std) + (1 | id),
     family = binomial("logit")),
  cc$frail_logratio, prior_binomial
)

fit_model(
  "cc_function_logratio",
  bf(mvbind(dinamometer_average_std, sft_cst_std, sft_elbow_flexion_std, sft_fut_std) ~
       time * (ilr_positive_to_positive_std + sex + age_std + imc_std) + (1 | id), family = gaussian()) +
    set_rescor(TRUE),
  cc$function_logratio, prior_for(fn_resp)
)
