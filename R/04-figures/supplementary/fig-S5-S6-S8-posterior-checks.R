# Posterior predictive checks
#   Figure S5: smoothed densities for the main models (FRAIL score and the four physical-function outcomes)
#   Figure S6: observed vs predicted counts at each FRAIL score (main, beta-binomial and ordinal models);
#              counts are also saved as a table (output/tables/diagnostics/ppc-frail-counts.csv)
#   Figure S8: smoothed densities for the imputation-based models (requires R/02-models/03-imputation-models.R)

source(here::here("R/functions.R"))
library(brms)
library(ggplot2)
library(patchwork)

out_dir <- here("manuscript/figures/supplementary/models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(1234)

# Figure S5 ---------------------------------------------------------------------------------------------

frail_fit <- readRDS(here("models/main_frail.rds"))
function_fit <- readRDS(here("models/main_function.rds"))
n_draws <- 1000

fig_s5 <- (pp_check(frail_fit, ndraws = n_draws) + ggtitle("FRAIL")) /
  (
    (pp_check(function_fit, resp = "dinamometeraveragestd", ndraws = n_draws) + ggtitle("Handgrip")) +
      (pp_check(function_fit, resp = "sftcststd", ndraws = n_draws) + ggtitle("CSTS")) +
      (pp_check(function_fit, resp = "sftelbowflexionstd", ndraws = n_draws) + ggtitle("Elbow flexion")) +
      (pp_check(function_fit, resp = "sftfutstd", ndraws = n_draws) + ggtitle("TUG"))
  ) +
  plot_layout(guides = "collect", heights = c(1/3, 2/3)) &
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "fig-S5.png"), plot = fig_s5, width = 850, height = 600,
       scale = 5, dpi = 450, units = "px")

# Figure S6 ---------------------------------------------------------------------------------------------

time_label <- c("t-1" = "Baseline", "t-2" = "6 months", "t-3" = "1 year")
n_draws_ppc <- 4000

get_counts <- function(name, label, levels, response) {
  fit <- readRDS(here("models", paste0(name, ".rds")))
  d <- as.data.table(fit$data)
  y <- d[[response]]
  if (is.factor(y)) y <- as.integer(y)
  yrep <- posterior_predict(fit, ndraws = n_draws_ppc)
  strata <- factor(time_label[as.character(d$time)], levels = time_label)
  ppc_counts(yrep, y, levels = levels, strata = strata)[, model := label][]
}

counts <- rbindlist(list(
  get_counts("main_frail", "Binomial (main)", 0:5, "frail_total"),
  get_counts("cc_frail_bb", "Beta-binomial", 0:5, "frail_total"),
  get_counts("cc_frail_ord", "Ordinal (0, 1, >=2)", 1:3, "frail_3cat")
))
counts[model == "Ordinal (0, 1, >=2)", score_label := c("0", "1", ">=2")[score]]
counts[model != "Ordinal (0, 1, >=2)", score_label := as.character(score)]

dir.create(here("output/tables/diagnostics"), recursive = TRUE, showWarnings = FALSE)
fwrite(counts, here("output/tables/diagnostics/ppc-frail-counts.csv"))

plot_data <- copy(counts)
plot_data[, stratum := factor(stratum, levels = c("Overall", time_label))]
plot_data[, score_label := factor(score_label, levels = c(as.character(0:5), ">=2"))]
plot_data[, model := factor(model, levels = c("Binomial (main)", "Beta-binomial", "Ordinal (0, 1, >=2)"))]

fig_s6 <- ggplot(plot_data, aes(x = score_label)) +
  geom_col(aes(y = observed), fill = "grey80", color = "grey40", width = 0.7) +
  geom_pointrange(aes(y = predicted, ymin = lower, ymax = upper), color = "#d62728", size = 0.3) +
  facet_grid(model ~ stratum, scales = "free", space = "free_x") +
  labs(x = "FRAIL score", y = "Number of observations",
       caption = "Bars: observed counts. Points and lines: posterior predictive mean and 95% interval.") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = NA))
ggsave(file.path(out_dir, "fig-S6.png"), plot = fig_s6, width = 3000, height = 2400, units = "px", dpi = 350)

# Figure S8 ---------------------------------------------------------------------------------------------

if (file.exists(here("models/mi_frail.rds")) && file.exists(here("models/mi_function.rds"))) {
  mi_frail <- readRDS(here("models/mi_frail.rds"))
  mi_function <- readRDS(here("models/mi_function.rds"))

  fig_s8 <- (pp_check(mi_frail, resp = "frailtotal", ndraws = n_draws) + ggtitle("FRAIL")) /
    (
      (pp_check(mi_function, resp = "dinamometeraveragestd", ndraws = n_draws) + ggtitle("Handgrip")) +
        (pp_check(mi_function, resp = "sftcststd", ndraws = n_draws) + ggtitle("CSTS")) +
        (pp_check(mi_function, resp = "sftelbowflexionstd", ndraws = n_draws) + ggtitle("Elbow flexion")) +
        (pp_check(mi_function, resp = "sftfutstd", ndraws = n_draws) + ggtitle("TUG"))
    ) / (
      (pp_check(mi_function, resp = "cd41positivecd62negativestd", ndraws = n_draws) + ggtitle("CD41+/CD62—")) +
        (pp_check(mi_function, resp = "cd41positivecd62positivestd", ndraws = n_draws) + ggtitle("CD41+/CD62+"))
    ) +
    plot_layout(guides = "collect", heights = c(1/3, 2/3, 1/3)) &
    theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "fig-S8.png"), plot = fig_s8, width = 850, height = 800,
         scale = 5, dpi = 450, units = "px")
} else {
  message("Imputation-based models not found: Figure S8 skipped.")
}
