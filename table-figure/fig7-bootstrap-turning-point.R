################################################################################
# Figure 8. Cluster-bootstrap distributions of the turning point h*
#            for the NB trip-count and ordered-logit trip-bin models
#
# Produces a two-panel figure:
#   (a) NB quadratic model  — total daily trip count
#   (b) Ordered logit (PO)  — trip-count bins
#
# Bootstrap: person-cluster resampling (preserves within-person dependence)
# B = 500 replications per model
#
# Depends on: 00-common-data.R  (provides trip_count, dat_psm)
################################################################################

library(ggplot2)
library(MASS)
library(ordinal)
library(dplyr)
library(patchwork)

set.seed(42)
B <- 500

################################################################################
# HELPERS
################################################################################

# Person-cluster bootstrap: resample person_ids with replacement
cluster_resample <- function(data, id_col = "person_id") {
  ids  <- unique(data[[id_col]])
  G    <- length(ids)
  by_id <- split(data, data[[id_col]])

  ids_b <- sample(ids, size = G, replace = TRUE)
  dat_list <- vector("list", G)
  for (k in seq_along(ids_b)) {
    d <- by_id[[as.character(ids_b[k])]]
    d$.boot_id <- k
    dat_list[[k]] <- d
  }
  dplyr::bind_rows(dat_list)
}

################################################################################
# PANEL (a): NB QUADRATIC — total trip count
################################################################################

# Point estimate
nb_model_q <- glm.nb(
  count ~ travel_dow + age + gender + education +
    telework_hour + I(telework_hour^2) +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  data = trip_count,
  weights = trip_count$day_weight
)

get_hstar_nb <- function(mod) {
  b1 <- coef(mod)["telework_hour"]
  b2 <- coef(mod)["I(telework_hour^2)"]
  if (is.na(b2) || b2 == 0) return(NA_real_)
  as.numeric(-b1 / (2 * b2))
}

hstar_nb_hat <- get_hstar_nb(nb_model_q)

# Bootstrap
hstar_nb_boot <- rep(NA_real_, B)
pb <- txtProgressBar(min = 0, max = B, style = 3)

for (b in seq_len(B)) {
  dat_b <- cluster_resample(trip_count)

  mod_b <- try(
    glm.nb(
      count ~ travel_dow + age + gender + education +
        telework_hour + I(telework_hour^2) +
        num_people + num_vehicles + res_type + income_broad +
        thrive_community_type + travel_date_season + year,
      data = dat_b,
      weights = dat_b$day_weight
    ),
    silent = TRUE
  )

  if (!inherits(mod_b, "try-error")) {
    hstar_nb_boot[b] <- get_hstar_nb(mod_b)
  }
  setTxtProgressBar(pb, b)
}
close(pb)

hstar_nb_boot <- hstar_nb_boot[is.finite(hstar_nb_boot)]
ci_nb <- quantile(hstar_nb_boot, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

cat(sprintf("\n--- NB count model ---\n"))
cat(sprintf("h* point estimate:  %.3f hrs\n", hstar_nb_hat))
cat(sprintf("95%% bootstrap CI:   [%.3f, %.3f]\n", ci_nb[1], ci_nb[3]))
cat(sprintf("Bootstrap median:   %.3f\n", ci_nb[2]))
cat(sprintf("Successful draws:   %d / %d\n\n", length(hstar_nb_boot), B))


################################################################################
# PANEL (b): ORDERED LOGIT (PO) — trip-count bins
################################################################################

# Standardize telework
trip_count_cz <- trip_count %>%
  mutate(
    telework_c  = telework_hour - mean(telework_hour, na.rm = TRUE),
    telework_z  = as.numeric(scale(telework_c)),
    telework_z2 = telework_z^2
  )

h_mean <- mean(trip_count_cz$telework_hour, na.rm = TRUE)
h_sd   <- sd(trip_count_cz$telework_hour,   na.rm = TRUE)

lvl <- levels(trip_count_cz$count_bin)

get_hstar_ord <- function(mod, h_mean, h_sd) {
  b <- coef(mod)
  if (!all(c("telework_z", "telework_z2") %in% names(b))) return(NA_real_)
  if (is.na(b["telework_z2"]) || b["telework_z2"] == 0)   return(NA_real_)
  z_star <- -b["telework_z"] / (2 * b["telework_z2"])
  as.numeric(h_mean + z_star * h_sd)
}

# Point estimate
ord_mod_z <- clm(
  count_bin ~ travel_dow + age + gender + education +
    telework_z + telework_z2 +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  data = trip_count_cz,
  weights = day_weight,
  link = "logit",
  Hess = TRUE
)

hstar_ord_hat <- get_hstar_ord(ord_mod_z, h_mean, h_sd)

# Bootstrap
hstar_ord_boot <- rep(NA_real_, B)
pb <- txtProgressBar(min = 0, max = B, style = 3)

for (b in seq_len(B)) {
  dat_b <- cluster_resample(trip_count_cz)
  dat_b$count_bin <- factor(dat_b$count_bin, levels = lvl, ordered = TRUE)

  # Recompute z and z2 within bootstrap sample
  dat_b <- dat_b %>%
    mutate(
      telework_c  = telework_hour - mean(telework_hour, na.rm = TRUE),
      telework_z  = as.numeric(scale(telework_c)),
      telework_z2 = telework_z^2
    )

  h_mean_b <- mean(dat_b$telework_hour, na.rm = TRUE)
  h_sd_b   <- sd(dat_b$telework_hour,   na.rm = TRUE)

  mod_b <- try(
    clm(
      count_bin ~ travel_dow + age + gender + education +
        telework_z + telework_z2 +
        num_people + num_vehicles + res_type + income_broad +
        thrive_community_type + travel_date_season + year,
      data = dat_b,
      weights = day_weight,
      link = "logit",
      Hess = TRUE
    ),
    silent = TRUE
  )

  if (!inherits(mod_b, "try-error")) {
    hstar_ord_boot[b] <- get_hstar_ord(mod_b, h_mean_b, h_sd_b)
  }
  setTxtProgressBar(pb, b)
}
close(pb)

hstar_ord_boot <- hstar_ord_boot[is.finite(hstar_ord_boot)]
ci_ord <- quantile(hstar_ord_boot, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

cat(sprintf("--- Ordered logit (bins) model ---\n"))
cat(sprintf("h* point estimate:  %.3f hrs\n", hstar_ord_hat))
cat(sprintf("95%% bootstrap CI:   [%.3f, %.3f]\n", ci_ord[1], ci_ord[3]))
cat(sprintf("Bootstrap median:   %.3f\n", ci_ord[2]))
cat(sprintf("Successful draws:   %d / %d\n\n", length(hstar_ord_boot), B))


################################################################################
# FIGURE: TWO-PANEL BOOTSTRAP HISTOGRAMS
################################################################################

# Build data frames for ggplot
df_nb <- data.frame(
  hstar = hstar_nb_boot,
  model = "NB trip-count model"
)

df_ord <- data.frame(
  hstar = hstar_ord_boot,
  model = "Ordered logit trip-bin model"
)

# Shared theme
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    panel.grid.minor = element_blank(),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 9)
  )

# Panel (a): NB count model
pa <- ggplot(df_nb, aes(x = hstar)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 35, fill = "#4E79A7", color = "white",
                 alpha = 0.7) +
  geom_density(linewidth = 0.6, color = "#2B5F8A") +
  geom_vline(xintercept = hstar_nb_hat,
             linewidth = 0.9, color = "black") +
  geom_vline(xintercept = ci_nb[c(1, 3)],
             linewidth = 0.7, linetype = "dashed", color = "#E15759") +
  annotate("text",
           x = hstar_nb_hat, y = Inf,
           label = sprintf("h* = %.2f", hstar_nb_hat),
           vjust = 2, hjust = -0.1,
           size = 3.5, fontface = "bold") +
  annotate("text",
           x = mean(c(ci_nb[1], ci_nb[3])), y = Inf,
           label = sprintf("95%% CI: [%.2f, %.2f]", ci_nb[1], ci_nb[3]),
           vjust = 3.5, size = 3, color = "#E15759") +
  labs(
    title    = "(a) Negative binomial model (total daily trip count)",
    x        = "Turning point h* (telework hours)",
    y        = "Density"
  ) +
  theme_pub

# Panel (b): Ordered logit bins model
pb_plot <- ggplot(df_ord, aes(x = hstar)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 35, fill = "#F28E2B", color = "white",
                 alpha = 0.7) +
  geom_density(linewidth = 0.6, color = "#C06D10") +
  geom_vline(xintercept = hstar_ord_hat,
             linewidth = 0.9, color = "black") +
  geom_vline(xintercept = ci_ord[c(1, 3)],
             linewidth = 0.7, linetype = "dashed", color = "#E15759") +
  annotate("text",
           x = hstar_ord_hat, y = Inf,
           label = sprintf("h* = %.2f", hstar_ord_hat),
           vjust = 2, hjust = -0.1,
           size = 3.5, fontface = "bold") +
  annotate("text",
           x = mean(c(ci_ord[1], ci_ord[3])), y = Inf,
           label = sprintf("95%% CI: [%.2f, %.2f]", ci_ord[1], ci_ord[3]),
           vjust = 3.5, size = 3, color = "#E15759") +
  labs(
    title    = "(b) Ordered logit model (trip-count bins)",
    x        = "Turning point h* (telework hours)",
    y        = "Density"
  ) +
  theme_pub

# Combine
fig8 <- pa / pb_plot +
  plot_annotation(
    title    = "Bootstrap distributions of the telework turning point (h*)",
    subtitle = "Person-cluster bootstrap, B = 500; solid line = point estimate, dashed lines = 95% percentile CI",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

fig8

ggsave("modeling/table-figure/fig8-bootstrap-turning-point.png",
       fig8, width = 8, height = 8, dpi = 300)

################################################################################
# SUMMARY TABLE (console output)
################################################################################

summary_table <- data.frame(
  Model         = c("NB trip-count", "Ordered logit trip-bin"),
  h_star        = round(c(hstar_nb_hat, hstar_ord_hat), 3),
  CI_lower      = round(c(ci_nb[1], ci_ord[1]), 3),
  CI_upper      = round(c(ci_nb[3], ci_ord[3]), 3),
  Median        = round(c(ci_nb[2], ci_ord[2]), 3),
  SE_boot       = round(c(sd(hstar_nb_boot), sd(hstar_ord_boot)), 3),
  N_success     = c(length(hstar_nb_boot), length(hstar_ord_boot)),
  stringsAsFactors = FALSE
)

cat("\n========================================\n")
cat("  BOOTSTRAP TURNING-POINT SUMMARY\n")
cat("========================================\n")
print(summary_table, row.names = FALSE)

# save data as all data right now as rdata
save.image("modeling/table-figure/fig8-bootstrap-turning-point-data.rdata")
