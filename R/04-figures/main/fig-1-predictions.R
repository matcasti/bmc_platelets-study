# Figure 1: model-based association of the platelet percentages with frailty and physical function
#
# Predictions from the main (complete-case) models, over the observed range of the platelet data.
# Prediction settings (as stated in the figure caption):
#   - one platelet percentage varies while the other is fixed at its mean (0, standardized units);
#     age and BMI are fixed at their means
#   - population-level predictions (random intercepts excluded)
#   - marginalized over time and sex, weighting each visit-by-sex combination by its observed proportion
#     among the visits used to fit each model
#   - frailty on the expected FRAIL score (0-5) scale; function outcomes in standardized units

source(here::here("R/functions.R"))
library(brms)
library(tidybayes)
library(ggplot2)

frail_fit <- readRDS(here("models/main_frail.rds"))
function_fit <- readRDS(here("models/main_function.rds"))

set.seed(1234)
n_draws <- 2000
x_step <- 0.1
negative <- "cd41_positive_cd62_negative_std"
positive <- "cd41_positive_cd62_positive_std"

predict_marginal <- function(fit, platelet, range_x) {
  d <- as.data.table(fit$data)
  w <- d[, .N, list(time, sex)][, w := N / sum(N)]
  other <- setdiff(c(negative, positive), platelet)

  grid <- CJ(x = seq(range_x[1], range_x[2], x_step), time = levels(d$time), sex = levels(d$sex))
  grid <- w[grid, on = c("time", "sex")][is.na(w), w := 0]
  newdata <- data.frame(time = factor(grid$time, levels = levels(d$time), ordered = TRUE),
                        sex = factor(grid$sex, levels = levels(d$sex)), age_std = 0, imc_std = 0)
  newdata[[platelet]] <- grid$x
  newdata[[other]] <- 0

  ep <- posterior_epred(fit, newdata = newdata, re_formula = NA, ndraws = n_draws)
  if (length(dim(ep)) == 2) ep <- array(ep, dim = c(dim(ep), 1), dimnames = list(NULL, NULL, "frailtotal"))
  out <- rbindlist(lapply(seq_len(dim(ep)[3]), function(k) {
    m <- ep[, , k]
    rbindlist(lapply(unique(grid$x), function(xv) {
      idx <- which(grid$x == xv)
      data.table(draw = seq_len(nrow(m)), z_score = xv, outcome = dimnames(ep)[[3]][k],
                 .epred = as.numeric(m[, idx, drop = FALSE] %*% grid$w[idx]) / sum(grid$w[idx]))
    }))
  }))
  out[, plaqueta := factor(if (platelet == negative) "CD41+CD62-" else "CD41+CD62+")]
  out[]
}

observed_range <- function(fit) range(c(as.data.table(fit$data)[[negative]], as.data.table(fit$data)[[positive]]))
range_of_data <- range(observed_range(frail_fit), observed_range(function_fit))

add_ribbons <- function(p) {
  for (w in c(0.95, 0.9, 0.8, 0.5)) {
    p <- p + stat_summary(geom = "ribbon", fun.data = median_qi, aes(fill = plaqueta), alpha = 0.2,
                          show.legend = FALSE, fun.args = list(.width = w))
  }
  p
}
style <- list(theme_classic(), scale_x_continuous(expand = c(0, 0)),
              scale_color_brewer(aesthetics = c("color", "fill"), type = "qual", palette = 6, direction = -1))

# Panel A: frailty -------------------------------------------------------------------------------------

epred_a <- rbind(predict_marginal(frail_fit, negative, range_of_data), predict_marginal(frail_fit, positive, range_of_data))
fig1a <- ggplot(epred_a, aes(z_score, .epred, group = plaqueta)) +
  facet_wrap(~plaqueta, nrow = 1) +
  stat_summary(geom = "line", fun = median, aes(color = plaqueta), linewidth = 1, show.legend = FALSE)
fig1a <- add_ribbons(fig1a) + style + labs(x = NULL, y = "Frailty score")

# Panel B: physical function ------------------------------------------------------------------------------

epred_b <- rbind(predict_marginal(function_fit, negative, range_of_data), predict_marginal(function_fit, positive, range_of_data))
epred_b[, outcome := factor(outcome, levels = c("dinamometeraveragestd", "sftcststd", "sftelbowflexionstd", "sftfutstd"),
                            labels = c("Dinamometer (kg)", "Chair Sit-to-Stand (reps)", "Elbow flexion (reps)", "TUG (seconds)"))]
fig1b <- ggplot(epred_b, aes(z_score, .epred, group = plaqueta)) +
  facet_grid(cols = vars(outcome), rows = vars(plaqueta), scales = "free_y") +
  stat_summary(geom = "line", fun = median, aes(color = plaqueta), linewidth = 1, show.legend = FALSE)
fig1b <- add_ribbons(fig1b) + style + labs(x = "Platelet Proportion\n(Standardized units)", y = "Standardized units")

dir.create(here("manuscript/figures/main"), recursive = TRUE, showWarnings = FALSE)
fig1 <- cowplot::plot_grid(fig1a, fig1b, align = "hv", axis = "lr", nrow = 2, rel_heights = c(0.3, 0.7), labels = "AUTO")
ggsave(here("manuscript/figures/main/fig-1.png"), plot = fig1, width = 500, height = 500, units = "px",
       scale = 7, dpi = 500, device = "png")
