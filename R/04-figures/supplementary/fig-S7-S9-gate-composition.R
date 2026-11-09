# Figure S7: gate composition by visit; Figure S9: percentage of CD41-positive events by assessment date
# (source tables: Tables S11 and S11b; see R/01-data/04-gate-composition.R)

source(here::here("R/functions.R"))
library(ggplot2)

load(here("data/platelets_d.RData"))
out_dir <- here("manuscript/figures/supplementary/gating")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- platelets_d[!is.na(cd41_positive_cd62_positive)]
d[, cd41_positive := cd41_positive_cd62_negative + cd41_positive_cd62_positive]
d[, cd62_in_cd41 := 100 * cd41_positive_cd62_positive / cd41_positive]
time_label <- c("t-1" = "Baseline", "t-2" = "6 months", "t-3" = "1 year")
d[, visit := factor(time_label[as.character(time)], levels = time_label)]

# Figure S7 ---------------------------------------------------------------------------------------------

plot_data <- melt(d, id.vars = c("id", "visit"), measure.vars = c("cd41_positive", "cd62_in_cd41"))
plot_data[, variable := factor(variable, levels = c("cd41_positive", "cd62_in_cd41"),
                               labels = c("CD41+ events (% of gated events)", "CD62+ events within the CD41+ gate (%)"))]

fig_s7 <- ggplot(plot_data, aes(visit, value)) +
  geom_line(aes(group = id), color = "grey70", linewidth = 0.3) +
  geom_boxplot(width = 0.35, outlier.shape = NA, fill = NA, linewidth = 0.5) +
  geom_point(position = position_jitter(width = 0.05, height = 0, seed = 1234), size = 1.2, alpha = 0.6, color = "#1f77b4") +
  facet_wrap(~ variable, scales = "free_y") +
  labs(x = NULL, y = "Percentage") +
  theme_classic(base_size = 12) +
  theme(strip.background = element_rect(fill = NA))
ggsave(file.path(out_dir, "fig-S7.png"), plot = fig_s7, width = 2800, height = 1400, units = "px", dpi = 350)

# Figure S9 ---------------------------------------------------------------------------------------------
# Assessment date is available as days since the start of data collection (the calendar origin is not distributed),
# so the axis shows months since the first assessment.

fig_s9 <- ggplot(d[!is.na(visit_day)], aes(visit_day / 30.44, cd41_positive, color = visit)) +
  geom_point(size = 1.3, alpha = 0.7) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Months since the start of data collection", y = "CD41+ events (% of cells gate)", color = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "fig-S9.png"), plot = fig_s9, width = 2000, height = 1200, units = "px", dpi = 300)
