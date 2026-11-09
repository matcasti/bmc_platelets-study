# Figure 2: platelet percentages across visits
#
# Model-based estimates from Gaussian mixed models of the observed values (models plt_neg and plt_pos),
# displayed with the observed values (lines join the samples of each participant).

source(here::here("R/functions.R"))
library(brms)
library(tidybayes)
library(ggplot2)

d <- prepare_analysis_data()
fits <- list("CD41+CD62-" = readRDS(here("models/plt_neg.rds")),
             "CD41+CD62+" = readRDS(here("models/plt_pos.rds")))
set.seed(1234)
time_labels <- c("Baseline", "6-Month", "1-Year")

epred_data <- rbindlist(lapply(names(fits), function(nm) {
  nd <- data.frame(time = factor(c("t-1", "t-2", "t-3"), levels = c("t-1", "t-2", "t-3"), ordered = TRUE))
  as.data.table(add_epred_draws(nd, fits[[nm]], re_formula = NA, ndraws = 2000))[
    , .category := nm][, time := factor(time, levels = c("t-1", "t-2", "t-3"), labels = time_labels)]
}))

dat <- melt.data.table(d[, c("time", "id", "cd41_positive_cd62_negative_std", "cd41_positive_cd62_positive_std")],
                       id.vars = c("time", "id"), value.name = "value", variable.name = ".category")
dat <- na.omit(dat)
dat[, .category := factor(.category, levels = c("cd41_positive_cd62_negative_std", "cd41_positive_cd62_positive_std"),
                          labels = c("CD41+CD62-", "CD41+CD62+"))]
dat[, time := factor(time, levels = c("t-1", "t-2", "t-3"), labels = time_labels)]

fig2 <- ggplot(epred_data, aes(time, .epred, group = .category)) +
  facet_grid(cols = vars(.category)) +
  stat_summary(geom = "line", fun = median, aes(color = .category), linewidth = 1, show.legend = FALSE)
for (w in c(0.95, 0.9, 0.8, 0.5)) {
  fig2 <- fig2 + stat_summary(geom = "ribbon", fun.data = median_qi, aes(fill = .category), alpha = 0.2,
                              show.legend = FALSE, fun.args = list(.width = w))
}
fig2 <- fig2 +
  geom_point(data = dat, aes(y = value, color = .category), cex = 0.3,
             position = position_jitter(width = 0.05, seed = 1234), show.legend = FALSE) +
  geom_line(data = dat, aes(y = value, color = .category, group = id), linewidth = 0.1, alpha = 0.2,
            position = position_jitter(width = 0.05, seed = 1234), show.legend = FALSE) +
  theme_classic() +
  scale_x_discrete(expand = c(0.2, 0)) +
  scale_color_brewer(aesthetics = c("color", "fill"), type = "qual", palette = 6, direction = -1) +
  labs(x = "Time", y = "Platelet proportion\n(Standardized units)")

dir.create(here("manuscript/figures/main"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("manuscript/figures/main/fig-2.png"), plot = fig2, width = 400, height = 300, units = "px",
       scale = 7, dpi = 500, device = "png")
