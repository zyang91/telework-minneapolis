################################################################################
# Figure 6. Covariate balance before and after matching and predicted
#            trip-bin probabilities in the matched sample
################################################################################


library(ggplot2)
library(patchwork)
library(ordinal)
library(tidyr)

################################################################################
# PANEL A: Covariate balance (standardized mean differences)
################################################################################

# Create telework_binary in main sample
trip_count <- trip_count %>%
  mutate(telework_binary = ifelse(telework_hour > 0, 1, 0))

# Covariates to check balance on
balance_vars <- c("age", "gender", "education", "num_people", "num_vehicles",
                  "income_broad", "thrive_community_type",
                  "travel_dow", "travel_date_season")

# Function to compute standardized mean difference (numeric or dummy-encoded)
compute_smd <- function(data, treat_var, covariates) {
  results <- list()
  for (v in covariates) {
    if (is.numeric(data[[v]])) {
      m1 <- mean(data[[v]][data[[treat_var]] == 1], na.rm = TRUE)
      m0 <- mean(data[[v]][data[[treat_var]] == 0], na.rm = TRUE)
      s1 <- sd(data[[v]][data[[treat_var]] == 1], na.rm = TRUE)
      s0 <- sd(data[[v]][data[[treat_var]] == 0], na.rm = TRUE)
      pooled_sd <- sqrt((s1^2 + s0^2) / 2)
      smd <- (m1 - m0) / pooled_sd
      results[[v]] <- data.frame(variable = v, smd = abs(smd))
    } else {
      levs <- unique(data[[v]])
      for (lv in levs) {
        m1 <- mean(data[[v]][data[[treat_var]] == 1] == lv, na.rm = TRUE)
        m0 <- mean(data[[v]][data[[treat_var]] == 0] == lv, na.rm = TRUE)
        pooled_sd <- sqrt((m1 * (1 - m1) + m0 * (1 - m0)) / 2)
        if (pooled_sd > 0) {
          smd <- (m1 - m0) / pooled_sd
        } else {
          smd <- 0
        }
        results[[paste0(v, ": ", lv)]] <-
          data.frame(variable = paste0(v, ": ", lv), smd = abs(smd))
      }
    }
  }
  do.call(rbind, results)
}

# Before matching
smd_before <- compute_smd(trip_count, "telework_binary", balance_vars) %>%
  mutate(sample = "Before matching")

# After matching
dat_psm <- dat_psm %>%
  mutate(telework_binary = as.numeric(as.character(telework_binary)),
         num_people   = as.numeric(num_people),
         num_vehicles = as.numeric(num_vehicles))

smd_after <- compute_smd(dat_psm, "telework_binary", balance_vars) %>%
  mutate(sample = "After matching")

smd_df <- bind_rows(smd_before, smd_after) %>%
  mutate(sample = factor(sample, levels = c("Before matching", "After matching")))

# Shorten long covariate labels for journal readability
label_map <- c(
  "education: Less than college degree" = "Edu: < College",
  "education: Bachelor's degree" = "Edu: Bachelor's",
  "education: Graduate/post-graduate degree" = "Edu: Graduate/post-grad",
  "education: Associate degree" = "Edu: Associate",
  "education: Missing" = "Edu: Missing",
  "thrive_community_type: Outside" = "Community: Outside",
  "thrive_community_type: Urban" = "Community: Urban",
  "thrive_community_type: Suburban" = "Community: Suburban",
  "thrive_community_type: Rural" = "Community: Rural",
  "income_broad: $100K+" = "Income: $100K+",
  "income_broad: $25-50K" = "Income: $25-50K",
  "income_broad: $50-75K" = "Income: $50-75K",
  "income_broad: <$25K" = "Income: <$25K",
  "income_broad: $75-100K" = "Income: $75-100K",
  "income_broad: Undisclosed" = "Income: Undisclosed",
  "travel_date_season: Summer" = "Season: Summer",
  "travel_date_season: Spring" = "Season: Spring",
  "travel_date_season: Fall" = "Season: Fall",
  "travel_date_season: Winter" = "Season: Winter",
  "travel_dow: Thursday" = "DOW: Thursday",
  "travel_dow: Tuesday" = "DOW: Tuesday",
  "travel_dow: Wednesday" = "DOW: Wednesday",
  "gender: Prefer not to answer" = "Gender: Prefer not to answer",
  "gender: Male" = "Gender: Male",
  "gender: Female" = "Gender: Female",
  "age: 35 to 44" = "Age: 35-44",
  "age: 55 to 64" = "Age: 55-64",
  "age: 18 to 24" = "Age: 18-24",
  "age: 16 to 17" = "Age: 16-17",
  "age: 65 to 74" = "Age: 65-74",
  "age: 25 to 34" = "Age: 25-34",
  "age: 75 or older" = "Age: 75+",
  "age: 45 to 54" = "Age: 45-54"
)

smd_df <- smd_df %>%
  mutate(variable = ifelse(variable %in% names(label_map),
                           label_map[variable], variable))

# Extract covariate group from label prefix and order by group then SMD
smd_df <- smd_df %>%
  mutate(cov_group = case_when(
    grepl("^Edu:", variable)       ~ "Education",
    grepl("^Income:", variable)    ~ "Income",
    grepl("^Community:", variable) ~ "Community",
    grepl("^Season:", variable)    ~ "Season",
    grepl("^DOW:", variable)       ~ "Day of week",
    grepl("^Gender:", variable)    ~ "Gender",
    grepl("^Age:", variable)       ~ "Age",
    TRUE                           ~ "Other"
  ))

# Compute max SMD per group (using "Before matching" for sorting groups)
group_order <- smd_df %>%
  filter(sample == "Before matching") %>%
  group_by(cov_group) %>%
  summarise(group_max_smd = max(smd), .groups = "drop") %>%
  arrange(group_max_smd) %>%
  mutate(group_rank = row_number())

# Use max SMD across both samples so all variables get a sort value
var_smd <- smd_df %>%
  group_by(variable, cov_group) %>%
  summarise(max_smd = max(smd), .groups = "drop") %>%
  left_join(group_order, by = "cov_group") %>%
  arrange(group_rank, max_smd)

var_order <- var_smd$variable

smd_df <- smd_df %>%
  mutate(variable = factor(variable, levels = var_order)) %>%
  left_join(group_order, by = "cov_group")

# Plot Panel A
p_a <- ggplot(smd_df, aes(x = smd, y = variable, color = sample, shape = sample)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Before matching" = "#E15759", "After matching" = "#4E79A7")) +
  scale_shape_manual(values = c("Before matching" = 1, "After matching" = 16)) +
  labs(x = "Absolute standardized mean difference",
       y = NULL,
       subtitle = "(a) Covariate balance before and after PSM",
       color = NULL, shape = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_text(size = 7))

################################################################################
# PANEL B: Predicted trip-bin probabilities in matched sample
################################################################################

# Standardize telework in PSM data
dat_psm <- dat_psm %>%
  mutate(
    telework_c  = telework_hour - mean(telework_hour, na.rm = TRUE),
    telework_z  = as.numeric(scale(telework_c)),
    telework_z2 = telework_z^2
  )

# Fit ordered logit on PSM data
ord_mod_psm <- clm(
  count_bin ~ travel_dow + age + gender + education +
    telework_z + telework_z2 +
    num_people + num_vehicles +
    income_broad + thrive_community_type + travel_date_season + year,
  data = dat_psm,
  link = "logit",
  Hess = TRUE
)

# Reference profile (PSM)
h_mean_psm <- mean(dat_psm$telework_hour, na.rm = TRUE)
h_sd_psm   <- sd(dat_psm$telework_hour, na.rm = TRUE)
c_mean_psm <- mean(dat_psm$telework_c, na.rm = TRUE)
c_sd_psm   <- sd(dat_psm$telework_c, na.rm = TRUE)

# --- Manual probability computation from clm coefficients ---
# Bypasses predict.clm entirely; builds the linear predictor by hand
# so we never call model.matrix on a single-row data frame.
clm_prob_grid <- function(mod, tw_hours, dat_psm, h_mean, c_mean, c_sd) {
  cf <- coef(mod)
  is_thr <- grepl("\\|", names(cf))
  theta  <- cf[is_thr]   # threshold parameters (K-1 values)
  beta   <- cf[!is_thr]  # slope parameters

  y_levels <- c("0", "1-2", "3-4", "5-6", "7+")
  K <- length(y_levels)

  # Build ONE full design matrix from the fitting data so we can grab
  # the column for the modal reference person.
  rhs <- reformulate(attr(terms(mod), "term.labels"))
  X_full <- model.matrix(rhs, dat_psm)
  X_full <- X_full[, names(beta), drop = FALSE]

  # Identify the "reference row": modal categorical + median continuous.
  # Pick the row closest to that profile (good enough for a reference person).
  # Faster: just compute eta for all rows and take the median eta as baseline,
  # then only vary telework terms.

  # Column indices for telework terms
  tw_cols  <- which(names(beta) %in% c("telework_z", "telework_z2"))
  # Non-telework part of eta (fixed at median across sample)
  eta_base_all <- X_full[, -tw_cols, drop = FALSE] %*% beta[-tw_cols]
  eta_base <- median(eta_base_all)

  results <- list()
  for (h in tw_hours) {
    z_val  <- ((h - h_mean) - c_mean) / c_sd
    z2_val <- z_val^2

    eta <- eta_base + beta["telework_z"] * z_val + beta["telework_z2"] * z2_val

    # Cumulative probabilities via logistic CDF
    CP <- plogis(theta - as.numeric(eta))

    P <- numeric(K)
    P[1] <- CP[1]
    for (j in 2:(K - 1)) P[j] <- CP[j] - CP[j - 1]
    P[K] <- 1 - CP[K - 1]

    results[[length(results) + 1]] <- c(h, P)
  }

  out <- as.data.frame(do.call(rbind, results))
  names(out) <- c("telework_hour", y_levels)
  out
}

tw_grid_psm <- seq(0, 8, by = 0.25)
pp_psm_df <- clm_prob_grid(ord_mod_psm, tw_grid_psm, dat_psm,
                            h_mean_psm, c_mean_psm, c_sd_psm)

plotdf_psm <- pp_psm_df %>%
  pivot_longer(-telework_hour, names_to = "bin", values_to = "prob") %>%
  mutate(bin = factor(bin, levels = c("0", "1-2", "3-4", "5-6", "7+")))

p_b <- ggplot(plotdf_psm, aes(x = telework_hour, y = prob, color = bin, linetype = bin)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = c("0" = "#4E79A7", "1-2" = "#F28E2B", "3-4" = "#59A14F",
               "5-6" = "#E15759", "7+" = "#B07AA1"),
    name = "Trip-count bin"
  ) +
  scale_linetype_manual(
    values = c("0" = "solid", "1-2" = "dashed", "3-4" = "dotdash",
               "5-6" = "longdash", "7+" = "twodash"),
    name = "Trip-count bin"
  ) +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Telework hours per day",
       y = "Predicted probability",
       subtitle = "(b) Predicted trip-bin probabilities in the PSM matched sample") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

fig6 <- p_a + p_b +
  plot_layout(widths = c(1, 1.2)) +
  plot_annotation(
    title = "Covariate balance and predicted trip-bin probabilities in the matched sample",
    theme = theme(plot.title = element_text(face = "bold", size = 11))
  )

fig6

ggsave("modeling/table-figure/fig6-psm-balance-bins.png",
       fig6, width = 13, height = 7, dpi = 300)

