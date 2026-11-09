# Conventional (frequentist) linear mixed models without any imputation (Table S9c)
#
# Same predictors as the main models (time by platelet percentages, sex, age and BMI; random
# intercept per participant), fitted to the visits with observed outcome, covariates and platelet
# percentages. Estimates are standardized coefficients with 95% Wald confidence intervals.

source(here::here("R/functions.R"))
library(lme4)

d <- prepare_analysis_data()
d[, `:=`(neg = cd41_positive_cd62_negative_std, pos = cd41_positive_cd62_positive_std)]

outcomes <- c("Handgrip strength" = "dinamometer_average_std", "Chair sit-to-stand" = "sft_cst_std",
              "Arm curl" = "sft_elbow_flexion_std", "Timed up-and-go" = "sft_fut_std")

res <- rbindlist(lapply(names(outcomes), function(label) {
  d[, y := get(outcomes[[label]])]
  dd <- d[complete.cases(d[, .(y, neg, pos, age_std, imc_std, sex)])]
  m <- lmer(y ~ time * (neg + pos + sex + age_std + imc_std) + (1 | id), data = dd, REML = TRUE)
  ci <- confint(m, method = "Wald")
  b <- fixef(m)
  rows <- c(neg = "CD41+CD62-", pos = "CD41+CD62+")
  data.table(Outcome = label, Term = rows, n_obs = nrow(dd), n_participants = uniqueN(dd$id),
             Estimate = sprintf("%.2f [%.2f, %.2f]", b[names(rows)], ci[names(rows), 1], ci[names(rows), 2]))
}))

dir.create(here("output/tables/models"), recursive = TRUE, showWarnings = FALSE)
fwrite(res, here("output/tables/models/s9c-mixed-models.csv"))
res
