# Observations per variable and visit (Table S5) and participants/observations entering each model
# (Table S7, S7b, S7c)

invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))  # keep non-ASCII characters (±, ⁺, ⁻) intact when run standalone
library(data.table)
library(here)

load(here("data/platelets_d.RData"))
out_dir <- here("output/tables/descriptives")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Table S5: number of non-missing observations per variable and visit ----------------------------------

s5_variables <- c(
  cd41_negative_cd62_negative = "CD41-CD62- proportion",
  cd41_negative_cd62_positive = "CD41-CD62+ proportion",
  cd41_positive_cd62_negative = "CD41+CD62- proportion",
  cd41_positive_cd62_positive = "CD41+CD62+ proportion",
  sex = "Sex", age = "Age",
  bp_systolic = "Systolic Blood Pressure", bp_diastolic = "Diastolic Blood Pressure",
  bp_pam = "Mean Arterial Pressure", bp_pp = "Pulse Pressure", fc = "Heart Rate",
  weight = "Weight", height = "Height", imc = "BMI", waist_hip_index = "Waist/Hip ratio",
  fat_total = "Body fat", muscle_total = "Muscle mass", body_water = "Body water", bone_mass = "Bone mass",
  frail_total = "FRAIL score", sft_cst = "CSTS", sft_elbow_flexion = "Elbow flexion test", sft_fut = "TUG",
  diabetes = "Diabetes", hypertension = "Hypertension", dyslipidemia = "Dyslipidemia",
  dinamometer_average = "Dinamometer (average)"
)

s5 <- melt(platelets_d[, c("time", "id", names(s5_variables)), with = FALSE],
           id.vars = c("time", "id"), value.factor = FALSE)[!is.na(value), .N, list(variable, time)] |>
  dcast(variable ~ time, value.var = "N")
s5 <- s5[match(names(s5_variables), variable)]
s5[, variable := unname(s5_variables[as.character(variable)])]
names(s5) <- c("Variable", "Baseline", "6-Month", "1-Year")
fwrite(s5, file.path(out_dir, "s5-observations.csv"))

# Table S7: participants and observations of each model family ---------------------------------------

d <- copy(platelets_d)
plt_ok <- !is.na(d$cd41_positive_cd62_negative) & !is.na(d$cd41_positive_cd62_positive)
covs_ok <- !is.na(d$sex) & !is.na(d$age) & !is.na(d$imc)

count <- function(ok, label, outcome) {
  data.table(
    Model = label, Outcome = outcome,
    `Observations (outcome and covariates available)` = sum(ok),
    `Participants (outcome and covariates available)` = uniqueN(d$id[ok]),
    `Observations (observed platelets)` = sum(ok & plt_ok),
    `Participants (observed platelets)` = uniqueN(d$id[ok & plt_ok]),
    `Baseline / 6-month / 1-year (observed platelets)` =
      paste(vapply(c("t-1", "t-2", "t-3"), function(t) sum(ok & plt_ok & d$time == t), numeric(1)), collapse = " / ")
  )
}

fn_vars <- c("dinamometer_average", "sft_cst", "sft_elbow_flexion", "sft_fut")
ok_function <- covs_ok & Reduce(`&`, lapply(fn_vars, function(v) !is.na(d[[v]])))
fwrite(rbind(
  count(covs_ok & !is.na(d$frail_total), "Frailty", "FRAIL score"),
  count(ok_function, "Physical function", "Handgrip, CSTS, arm curl and TUG (joint model)")
), file.path(out_dir, "s7-models.csv"))

# Paired platelet data
wide <- dcast(d, id ~ time, value.var = "cd41_positive_cd62_negative")
fwrite(data.table(
  Quantity = c("Participants", "Visits attended", "Visits with platelet data",
               "Participants with platelet data at baseline and 1 year",
               "Participants with platelet data at all three visits",
               "Participants with at least one platelet measurement"),
  N = c(uniqueN(d$id), nrow(d), sum(plt_ok),
        sum(!is.na(wide$`t-1`) & !is.na(wide$`t-3`)),
        sum(!is.na(wide$`t-1`) & !is.na(wide$`t-2`) & !is.na(wide$`t-3`)),
        uniqueN(d$id[plt_ok]))
), file.path(out_dir, "s7-paired.csv"))

# Completeness of the additional confounders
fwrite(data.table(
  Variable = c("Diabetes", "Statin use"),
  `Non-missing visits` = c(sum(!is.na(d$diabetes)), sum(!is.na(d$statins))),
  `Participants with the condition/drug at any visit` = c(d[diabetes == "Yes", uniqueN(id)], d[statins == "Yes", uniqueN(id)]),
  `Participants whose status changed between visits` = c(d[, uniqueN(diabetes) > 1, id][, sum(V1)],
                                                         d[, uniqueN(statins) > 1, id][, sum(V1)])
), file.path(out_dir, "s7-confounders.csv"))
