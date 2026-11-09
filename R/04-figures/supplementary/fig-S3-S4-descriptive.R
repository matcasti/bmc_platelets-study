# Figure S3: platelet subtypes, absolute platelet counts and chronological age
# Figure S4: patient-level measurement plot (visits attended by each participant)

source(here::here("R/functions.R"))
library(ggplot2)

load(here("data/platelets_d.RData"))
out_dir <- here("manuscript/figures/supplementary/descriptive")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Figure S3 ------------------------------------------------------------------------------------------

plot_data <- age_platelet_correlations(platelets_d)$data

fig_a <- ggplot(plot_data, aes(Age, Proportion, color = `Platelet subtype`, fill = `Platelet subtype`)) +
  facet_wrap(~ `Platelet subtype`, nrow = 2, scales = "free_y") +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_brewer(aesthetics = c("color", "fill"), type = "qual", palette = 2) +
  scale_y_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  labs(caption = "Relative, subtype-wise, platelet count") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

fig_b <- ggplot(plot_data, aes(Age, Platelets)) +
  geom_point(color = "darkblue") +
  geom_smooth(method = "lm", fill = "darkblue", color = "darkblue") +
  labs(y = "Thousands per µL", caption = "Absolute total platelet count") +
  theme_classic(base_size = 12)

fig_s3 <- cowplot::plot_grid(fig_a, fig_b, align = "hv", axis = "lr", nrow = 2,
                             rel_heights = c(0.6, 0.4), labels = "AUTO")
ggsave(file.path(out_dir, "fig-S3.png"), plot = fig_s3, width = 400, height = 500, units = "px",
       scale = 7, dpi = 500, device = "png")

# Figure S4 ------------------------------------------------------------------------------------------

d <- copy(platelets_d)
d[, id := as.factor(id)]
nob_per_id <- d[, list(nobs = factor(.N, levels = 1:3)), id]
d <- nob_per_id[d, on = "id"]

fig_s4 <- ggplot(data = d, mapping = aes(y = time, x = reorder(id, as.numeric(nobs)))) +
  geom_line(mapping = aes(col = nobs, group = id), show.legend = TRUE) +
  geom_point(mapping = aes(col = nobs, shape = nobs)) +
  scale_y_discrete(expand = c(0, 0.05)) +
  scale_x_discrete(expand = c(0, 2), labels = NULL) +
  scale_color_brewer(palette = 2, type = "qual") +
  theme_classic() +
  labs(x = "Record ID", y = "Measurement time", color = "Number of measurements", shape = "Number of measurements") +
  theme(panel.grid.major.x = element_line(colour = "#0002"), legend.position = "top")

ggsave(file.path(out_dir, "fig-S4.png"), plot = fig_s4, width = 850, height = 600,
       scale = 5, dpi = 450, units = "px")
