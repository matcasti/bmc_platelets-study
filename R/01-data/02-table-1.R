# Table 1: sample description by sex (first record of each participant)

library(data.table)
library(here)
library(gtsummary)

load(here("data/platelets_d.RData"))

tbl_variables <- c("sex", "age", "diabetes", "hypertension", "dyslipidemia",
                   "angiotensin_receptor_blockers", "calcium_channel_blockers",
                   "statins", "ace_inhibitor", "hypoglycemics", "beta_blockers",
                   "bp_systolic", "bp_diastolic", "bp_pam", "bp_pp",
                   "fc", "weight", "height", "imc", "fat_total",
                   "muscle_total", "muscle_apendicular", "muscle_squeletal_index", "body_water", "bone_mass",
                   "dinamometer_average", "dinamometer_left_hand", "dinamometer_right_hand", "frail_total",
                   "sft_cst", "sft_elbow_flexion", "sft_fut",
                   "cd41_negative_cd62_negative", "cd41_negative_cd62_positive",
                   "cd41_positive_cd62_negative", "cd41_positive_cd62_positive")

tbl_data <- platelets_d[, .SD[1], by = id][, .SD, .SDcols = tbl_variables]

tbl_data[, sex := `levels<-`(sex, c("Male", "Female"))]

# FRAIL score ranges from 0 to 5: keep empty levels so every score is displayed
tbl_data[, frail_total := factor(frail_total, levels = 0:5)]

# The FRAIL score is displayed by level, but the male-female difference is the
# standardized mean difference of the (numeric) score
smd_ordinal <- function(data, variable, by, ...) {
  res <- smd::smd(x = as.numeric(as.character(data[[variable]])), g = data[[by]],
                  na.rm = TRUE, std.error = TRUE)
  data.frame(estimate = res$estimate,
             conf.low = res$estimate - 1.96 * res$std.error,
             conf.high = res$estimate + 1.96 * res$std.error)
}

tbl_1 <- tbl_summary(
  tbl_data,
  by = sex,
  missing = "no",
  statistic = all_continuous() ~ "{mean} ± {sd}",
  digits = list(all_continuous() ~ 1, all_dichotomous() ~ 0),
  label = list(
    age ~ "Age (years old)",
    diabetes ~ "Diabetes (n)",
    hypertension ~ "Hypertension (n)",
    dyslipidemia ~ "Dyslipidemia (n)",
    angiotensin_receptor_blockers ~ "Angiotensin receptor blockers (n)",
    calcium_channel_blockers ~ "Calcium channel blockers (n)",
    statins ~ "Statins (n)",
    ace_inhibitor ~ "ACE inhibitor (n)",
    hypoglycemics ~ "Hypoglycemic drugs (n)",
    beta_blockers ~ "Beta blockers (n)",
    bp_systolic ~ "Systolic blood pressure (mmHg)",
    bp_diastolic ~ "Diastolic blood pressure (mmHg)",
    bp_pam ~ "Mean arterial pressure (mmHg)",
    bp_pp ~ "Pulse pressure (mmHg)",
    fc ~ "Heart rate (bpm)",
    weight ~ "Weight (kg)",
    height ~ "Height (cm)",
    imc ~ "BMI (kg/m²)",
    fat_total ~ "Total fat mass (%)",
    muscle_total ~ "Total muscle mass (kg)",
    muscle_apendicular ~ "Appendicular muscle mass (kg)",
    muscle_squeletal_index ~ "Skeletal muscle index (kg/m²)",
    body_water ~ "Body water (%)",
    bone_mass ~ "Bone mass (kg)",
    dinamometer_average ~ "HGS - Average (kg)",
    dinamometer_left_hand ~ "HGS - Left hand (kg)",
    dinamometer_right_hand ~ "HGS - Right hand (kg)",
    frail_total ~ "Total FRAIL score",
    sft_cst ~ "Chair sit-to-stand (reps)",
    sft_elbow_flexion ~ "Elbow flexion test (reps)",
    sft_fut ~ "Timed-up and Go (seconds)",
    cd41_negative_cd62_negative ~ "CD41⁻CD62⁻ (%)",
    cd41_negative_cd62_positive ~ "CD41⁻CD62⁺ (%)",
    cd41_positive_cd62_negative ~ "CD41⁺CD62⁻ (%)",
    cd41_positive_cd62_positive ~ "CD41⁺CD62⁺ (%)"
  )
) |>
  add_overall() |>
  add_difference(
    test = list(all_continuous() ~ "smd", frail_total ~ smd_ordinal),
    include = !c(diabetes, hypertension, dyslipidemia, angiotensin_receptor_blockers,
                 calcium_channel_blockers, statins, ace_inhibitor, hypoglycemics, beta_blockers)
  ) |>
  remove_footnote_header() |>
  remove_abbreviation()

dir.create(here("output/tables/descriptives"), recursive = TRUE, showWarnings = FALSE)
saveRDS(tbl_1, file = here("output/tables/descriptives/tbl-1.rds"))
