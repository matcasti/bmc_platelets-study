# Prepare the analysis dataset (data/platelets_d.RData) from the study database
#
# Input : data/raw/data-raw.csv  (REDCap export; contains direct identifiers and is NOT distributed)
# Output: data/platelets_d.RData (de-identified: no names, national IDs, contact data, free-text medication
#         or calendar dates; visit dates are kept only as days since the first assessment and as month)
#
# The analysis dataset is distributed with the repository, so this script only needs to be run by the
# data custodians. Every other script starts from data/platelets_d.RData.

library(data.table)
library(here)

raw_file <- here("data/raw/data-raw.csv")
if (!file.exists(raw_file)) stop("data/raw/data-raw.csv not found: the raw database is not distributed.")

# Variables read from the database ------------------------------------------------------

phases <- 1:4
quadrants <- c("cd41_negative_cd62_positive", "cd41_positive_cd62_positive",
               "cd41_positive_cd62_negative", "cd41_negative_cd62_negative")

selected <- c(
  "record_id", "redcap_event_name", "sex", "age", "phase_date", "disease_physical", "drugs",
  "bp_systolic", "bp_diastolic", "bp_pam", "bp_pp", "fc", "weight", "height", "imc", "waist_hip_index",
  "fat_total", "muscle_total", "muscle_left_arm", "muscle_right_arm", "muscle_right_leg", "muscle_left_leg",
  "muscle_torso", "body_water", "bone_mass", "dinamometer_left_hand", "dinamometer_right_hand",
  "frail_total", "sft_cst", "sft_elbow_flexion", "sft_fut", "plt_10_9_l",
  as.vector(outer(c(quadrants, "aba_fcs_singlets", "aba_fcs_singlets_cells"), phases, paste, sep = "_fase_"))
)

data_raw <- fread(raw_file, select = selected)
data_raw[, id := rleid(record_id)]
data_raw[, record_id := NULL]

# Flow cytometry variables are stored in the row of the corresponding phase ---------------------

fcs_pattern <- "^cd41_|^aba_fcs_"
fcs_raw <- melt(data_raw, id.vars = c("id", "redcap_event_name"), measure.vars = patterns(fcs_pattern))
fcs_raw[, variable := gsub("_fase_[1-4]", "", variable)]
fcs_raw[, value := as.numeric(gsub(",", ".", value))]
fcs_raw <- na.omit(fcs_raw)
fcs_wide <- dcast(fcs_raw, formula = id + redcap_event_name ~ variable, value.var = "value")

data_raw[, grep("_fase_[1-4]$", names(data_raw), value = TRUE) := NULL]
data_merged <- fcs_wide[data_raw, on = c("id", "redcap_event_name")]
data_merged <- data_merged[redcap_event_name != "initial_data_arm_1"]

# Visit timing (kept only as days since the first assessment and as year-month) ------------------

data_merged[, phase_date := as.Date(phase_date)]
# one baseline date precedes the start of the study (data-entry error): treated as missing
data_merged[phase_date < as.Date("2022-07-01"), phase_date := NA]
data_merged[, visit_day := as.integer(phase_date - min(phase_date, na.rm = TRUE))]
data_merged[, visit_month := format(phase_date, "%Y-%m")]
data_merged[, phase_date := NULL]

setnames(data_merged, c("aba_fcs_singlets", "aba_fcs_singlets_cells"), c("singlets_pct", "cells_pct"))

# Factors ----------------------------------------------------------------------------------------

data_merged[, sex := factor(sex, levels = c("M", "F"))]

phase_labels <- c("fase_1_arm_1", "fase_2_arm_1", "fase_3_arm_1", "fase_4_arm_1")
data_merged[, year := factor(redcap_event_name, levels = phase_labels, labels = c("Year_1", "Year_1", "Year_2", "Year_2"))]
data_merged[, semester := factor(redcap_event_name, levels = phase_labels,
                                 labels = c("Semester_1", "Semester_2", "Semester_1", "Semester_2"))]
# t-0 is the first phase of the database, which is not part of the analyses (see below)
data_merged[, time := factor(redcap_event_name, levels = phase_labels,
                             labels = c("t-0", "t-1", "t-2", "t-3"), ordered = TRUE)]
data_merged[, redcap_event_name := NULL]

# Comorbidities and medication (derived from free text, which is not distributed) ---------------------

data_merged[, `:=`(diabetes = factor("No"), hypertension = factor("No"), dyslipidemia = factor("No"))]
data_merged[tolower(disease_physical) %like% "dm|insulin|diabet", diabetes := "Yes"]
data_merged[tolower(disease_physical) %like% "hta|hiperten|arterial", hypertension := "Yes"]
# hypercholesterolemia, dyslipidemia, triglycerides (spelling variants included)
data_merged[tolower(disease_physical) %like% "dislip|coleste|trigli|hiperlip", dyslipidemia := "Yes"]
data_merged[, disease_physical := NULL]

# a condition recorded at any visit is assigned to all visits of the participant
for (v in c("diabetes", "hypertension", "dyslipidemia")) {
  ids <- data_merged[, any(get(v) == "Yes"), id][V1 == TRUE, id]
  data_merged[id %in% ids, (v) := "Yes"]
}

data_merged[, drugs := tolower(drugs)]
drug_patterns <- c(angiotensin_receptor_blockers = "lozart|losart",
                   calcium_channel_blockers = "lodipino|nifedipino",
                   statins = "atorvasta|rubastat|statin",
                   ace_inhibitor = "enalapril",
                   hypoglycemics = "metformina",
                   beta_blockers = "carvedilol|propanolol|bisoprolol|atenolol")
for (v in names(drug_patterns)) {
  data_merged[, (v) := factor(fifelse(drugs %like% drug_patterns[[v]], "Yes", "No"))]
}
data_merged[, drugs := NULL]
for (v in c("diabetes", "hypertension", "dyslipidemia")) set(data_merged, j = v, value = factor(data_merged[[v]]))

# Cleaning of recorded values ----------------------------------------------------------------------------

data_merged[muscle_left_arm %like% "[0-9]\\,[0-9]|[0-9]\\/[0-9]", muscle_left_arm := gsub("\\,|\\/", "\\.", muscle_left_arm)]
data_merged[muscle_left_arm %like% "\\.$", muscle_left_arm := gsub("\\.$", "", muscle_left_arm)]
data_merged[muscle_left_arm == "", muscle_left_arm := NA_character_]
data_merged[, muscle_left_arm := as.numeric(muscle_left_arm)]
data_merged[muscle_left_arm > 10, muscle_left_arm := muscle_left_arm / 10]
data_merged[muscle_right_arm > 10, muscle_right_arm := muscle_right_arm / 10]

data_merged[, sft_fut := trimws(sft_fut)]
data_merged[sft_fut == "" | tolower(sft_fut) == "n/a", sft_fut := NA_character_]
data_merged[, sft_fut := gsub(" segundos", "", sft_fut)]
data_merged[, sft_fut := as.numeric(gsub("\\,|\\:", "\\.", sft_fut))]
data_merged[sft_fut > 20, sft_fut := sft_fut / 10]

# sex missing in some visits: filled from the other visits of the same participant
data_merged[id == 27, sex := "F"]
data_merged[id %in% c(91, 104), sex := "M"]

# Analysis visits: phases 2-4 (t-1, t-2, t-3), participants with at least one of them ------------------------
platelets_d <- data_merged[time != "t-0"]

platelets_d[height == 60, height := 160]
platelets_d[height < 100, height := height * 100]
platelets_d[, imc := weight / (height / 100)^2]

platelets_d[body_water > 100, body_water := body_water / 10]
platelets_d[body_water < 30 & id == 24, body_water := 46.05]
platelets_d[body_water < 30 & id == 41, body_water := 39.95]

platelets_d[muscle_total > 100, muscle_total := muscle_total / 100]
platelets_d[bone_mass > 5, bone_mass := bone_mass / 10]

platelets_d[, muscle_apendicular := muscle_left_arm + muscle_right_arm + muscle_left_leg + muscle_right_leg]
platelets_d[, muscle_squeletal_index := muscle_apendicular / (height / 100)^2]
platelets_d[, dinamometer_average := (dinamometer_left_hand + dinamometer_right_hand) / 2]

platelets_d <- droplevels(platelets_d)

save(platelets_d, file = here("data/platelets_d.RData"))
