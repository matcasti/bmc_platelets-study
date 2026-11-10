# Empirical bootstrap contrasts of the platelet percentages between visits (Table S4)
#
# Paired within-participant differences between visits, bootstrap of the mean (10,000 resamples).

invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))  # keep non-ASCII characters (±, ⁺, ⁻) intact when run standalone
library(data.table)
library(here)

load(here("data/platelets_d.RData"))
set.seed(1234)

long_platelets <- melt.data.table(
  platelets_d, id.vars = c("id", "time"),
  measure.vars = c("cd41_positive_cd62_negative", "cd41_positive_cd62_positive"),
  variable.name = "Phenotype")

wide_platelets <- dcast(long_platelets, formula = id + Phenotype ~ time)

paired_differences <- wide_platelets[, list(
  "Baseline to 6-Months" = `t-2` - `t-1`,
  "Baseline to 1-Year" = `t-3` - `t-1`,
  "6-Months to 1-Year" = `t-3` - `t-2`
), keyby = list(id, Phenotype)]

paired_long <- melt.data.table(paired_differences, id.vars = c("id", "Phenotype"), measure.vars = 3:5)

nboot <- 10000
boot <- paired_long[, vapply(1:nboot, function(x) {
  value <- value[!is.na(value)]
  mean(sample(value, replace = TRUE))
}, 1), list(Phenotype, variable)]

results <- boot[, list(Estimate = mean(V1), SE = sd(V1),
                       Q2.5 = quantile(V1, 0.025), Q50 = median(V1), Q97.5 = quantile(V1, 0.975)),
                list(Phenotype, variable)]
results[, Phenotype := factor(Phenotype,
                              levels = c("cd41_positive_cd62_negative", "cd41_positive_cd62_positive"),
                              labels = c("CD41+/CD62—", "CD41+/CD62+"))]
setnames(results, "variable", "Comparison")

dir.create(here("output/tables/gate"), recursive = TRUE, showWarnings = FALSE)
fwrite(results, here("output/tables/gate/s4-platelet-changes.csv"))
