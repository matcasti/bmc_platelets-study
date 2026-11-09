# Composition of the cells gate by visit and by date of assessment (Table S11, S11b, S11c)
#
# The percentage of CD41-positive events (of the cells gate) and of CD62-positive events within the
# CD41-positive gate, by visit, with paired within-participant changes (bootstrap), and the
# relation with the date of assessment (the only available proxy of the assay batch).
# The analysis dataset stores the assessment date as days since the first assessment and as year-month.

library(data.table)
library(here)

load(here("data/platelets_d.RData"))
set.seed(1234)

out_dir <- here("output/tables/gate")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- platelets_d[!is.na(cd41_positive_cd62_positive)]
d[, cd41_positive := cd41_positive_cd62_negative + cd41_positive_cd62_positive]
d[, cd62_in_cd41 := 100 * cd41_positive_cd62_positive / cd41_positive]

time_label <- c("t-1" = "Baseline", "t-2" = "6 months", "t-3" = "1 year")
d[, visit := factor(time_label[as.character(time)], levels = time_label)]

# Table S11: median [IQR] by visit and paired differences ----------------------------------------------

measures <- c(cd41_positive_cd62_negative = "CD41+CD62- (% of gated events)",
              cd41_positive_cd62_positive = "CD41+CD62+ (% of gated events)",
              cd41_negative_cd62_negative = "CD41-CD62- (% of gated events)",
              cd41_negative_cd62_positive = "CD41-CD62+ (% of gated events)",
              cd41_positive = "CD41+ events (% of gated events)",
              cd62_in_cd41 = "CD62+ events within the CD41+ gate (%)")

q <- function(x) sprintf("%.1f [%.1f, %.1f]", median(x), quantile(x, 0.25), quantile(x, 0.75))

paired_change <- function(v, from, to) {
  w <- dcast(d, id ~ time, value.var = v)
  diffs <- (w[[to]] - w[[from]])[!is.na(w[[to]]) & !is.na(w[[from]])]
  boot <- replicate(10000, mean(sample(diffs, replace = TRUE)))
  sprintf("%.1f [%.1f, %.1f]", mean(diffs), quantile(boot, 0.025), quantile(boot, 0.975))
}
n_paired <- function(from, to) {
  w <- dcast(d, id ~ time, value.var = "cd41_positive")
  sum(!is.na(w[[from]]) & !is.na(w[[to]]))
}

s11 <- rbindlist(lapply(names(measures), function(v) {
  data.table(Measure = measures[[v]],
             Baseline = q(d[time == "t-1"][[v]]), `6 months` = q(d[time == "t-2"][[v]]), `1 year` = q(d[time == "t-3"][[v]]),
             `Paired change, baseline to 1 year` = paired_change(v, "t-1", "t-3"),
             `Paired change, baseline to 6 months` = paired_change(v, "t-1", "t-2"))
}))
s11 <- rbind(data.table(Measure = "Samples (n)", Baseline = as.character(d[time == "t-1", .N]),
                        `6 months` = as.character(d[time == "t-2", .N]), `1 year` = as.character(d[time == "t-3", .N]),
                        `Paired change, baseline to 1 year` = as.character(n_paired("t-1", "t-3")),
                        `Paired change, baseline to 6 months` = as.character(n_paired("t-1", "t-2"))), s11)
fwrite(s11, file.path(out_dir, "s11-by-visit.csv"))

fwrite(d[, list(`Samples with <50% CD41+ events` = sum(cd41_positive < 50),
                `Samples with <80% CD41+ events` = sum(cd41_positive < 80), n = .N), visit],
       file.path(out_dir, "s11-low-cd41.csv"))

# Table S11b/c: gate summaries and CD41+ share by assessment date ----------------------------------------------

month_date <- as.Date(paste0(d$visit_month, "-15"))
d[, month_date := month_date]

by_visit <- d[, list(
  Samples = .N,
  `Samples with a visit date` = sum(!is.na(visit_month)),
  `Median visit month` = format(median(month_date, na.rm = TRUE), "%B %Y"),
  `Range of visit months` = paste(format(min(month_date, na.rm = TRUE), "%b %Y"), format(max(month_date, na.rm = TRUE), "%b %Y"), sep = " - "),
  `Singlets (% of events)` = q(singlets_pct),
  `Cells (% of singlets)` = q(cells_pct),
  `CD41+ events (% of cells gate)` = q(cd41_positive),
  `Spearman rho, CD41+ share vs visit date` = {
    ok <- !is.na(visit_day)
    r <- suppressWarnings(cor.test(visit_day[ok], cd41_positive[ok], method = "spearman", exact = FALSE))
    sprintf("%.2f (p %s)", r$estimate, if (r$p.value < 0.001) "< 0.001" else sprintf("= %.3f", r$p.value))
  }), visit]
setnames(by_visit, "visit", "Visit")
fwrite(by_visit, file.path(out_dir, "s11-dates.csv"))

fwrite(d[!is.na(visit_month), list(Samples = .N, `CD41+ events (% of cells gate), median [IQR]` = q(cd41_positive)),
         list(Visit = visit, Month = visit_month)][order(Visit, Month)],
       file.path(out_dir, "s11-months.csv"))
