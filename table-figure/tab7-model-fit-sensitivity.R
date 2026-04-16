################################################################################
# Table 7. Model fit and sensitivity checks for trip-bin and count models
################################################################################


library(MASS)
library(sandwich)
library(lmtest)
library(survey)
library(ordinal)
library(gt)

################################################################################
# PANEL A: COUNT MODEL FIT — Poisson vs NB (linear) vs NB (quadratic)
################################################################################

# --- Poisson baseline ---
poisson_model <- glm(
  count ~ travel_dow + age + gender + education +
    telework_hour + num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  data = trip_count,
  weights = trip_count$day_weight,
  family = poisson(link = "log")
)

overdispersion_pear <- sum(residuals(poisson_model, type = "pearson")^2) /
  poisson_model$df.residual

# --- NB linear ---
nb_linear <- glm.nb(
  count ~ travel_dow + age + gender + education +
    telework_hour + num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  data = trip_count,
  weights = trip_count$day_weight
)

# --- NB quadratic ---
nb_quad <- glm.nb(
  count ~ travel_dow + age + gender + education +
    telework_hour + I(telework_hour^2) +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  data = trip_count,
  weights = trip_count$day_weight
)

# Turning point
b1 <- coef(nb_quad)["telework_hour"]
b2 <- coef(nb_quad)["I(telework_hour^2)"]
tp_nb <- -b1 / (2 * b2)

# Weighted observed vs predicted (calibration MAE for NB quad)
y_all   <- trip_count$count
w_all   <- trip_count$day_weight
theta_q <- nb_quad$theta
mu_q    <- fitted(nb_quad)
kvals   <- 0:10
obs_w   <- sapply(kvals, function(k) sum(w_all[y_all == k], na.rm = TRUE))
pred_w  <- sapply(kvals, function(k)
  sum(w_all * dnbinom(k, size = theta_q, mu = mu_q), na.rm = TRUE))
tot_w   <- sum(w_all, na.rm = TRUE)
mae_nb  <- mean(abs(obs_w / tot_w - pred_w / tot_w))

# Same for NB linear
mu_lin    <- fitted(nb_linear)
theta_lin <- nb_linear$theta
pred_w_lin <- sapply(kvals, function(k)
  sum(w_all * dnbinom(k, size = theta_lin, mu = mu_lin), na.rm = TRUE))
mae_nb_lin <- mean(abs(obs_w / tot_w - pred_w_lin / tot_w))

count_fit <- data.frame(
  Model = c("Poisson", "NB (linear telework)", "NB (quadratic telework)"),
  AIC   = c(AIC(poisson_model), AIC(nb_linear), AIC(nb_quad)),
  logLik = c(as.numeric(logLik(poisson_model)),
             as.numeric(logLik(nb_linear)),
             as.numeric(logLik(nb_quad))),
  N_params = c(length(coef(poisson_model)),
               length(coef(nb_linear)) + 1,
               length(coef(nb_quad)) + 1),
  theta = c(NA, nb_linear$theta, nb_quad$theta),
  Overdispersion = c(round(overdispersion_pear, 2), NA, NA),
  Calibration_MAE = c(NA, round(mae_nb_lin, 5), round(mae_nb, 5)),
  Turning_point = c(NA, NA, round(tp_nb, 2)),
  stringsAsFactors = FALSE
)

################################################################################
# PANEL B: ORDINAL BIN MODEL FIT — PO (svyolr) vs PPOM (clm nominal)
################################################################################

trip_count_cz <- trip_count %>%
  mutate(
    telework_c  = telework_hour - mean(telework_hour, na.rm = TRUE),
    telework_z  = as.numeric(scale(telework_c)),
    telework_z2 = telework_z^2
  )

# PO model (survey-weighted)
svy_design <- svydesign(ids = ~person_id, weights = ~day_weight, data = trip_count_cz)

ord_po <- svyolr(
  count_bin ~ travel_dow + age + gender + education +
    telework_z + telework_z2 +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  design = svy_design
)

# PPOM model (partial proportional odds, clm with nominal telework)
m_ppom <- clm(
  count_bin ~ travel_dow + age + gender + education +
    telework_z + telework_z2 +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  nominal = ~ telework_z,
  data = trip_count_cz,
  weights = day_weight,
  link = "logit",
  Hess = TRUE
)

# Weighted log scores for PO vs PPOM
P_po <- predict(ord_po, newdata = trip_count_cz, type = "probs")
wt_norm <- trip_count_cz$day_weight / sum(trip_count_cz$day_weight, na.rm = TRUE)

# Match observed bin to column index
bin_levels <- c("0", "1-2", "3-4", "5-6", "7+")
obs_idx <- match(as.character(trip_count_cz$count_bin), bin_levels)

p_obs_po <- sapply(seq_len(nrow(P_po)), function(i) P_po[i, obs_idx[i]])
p_obs_po <- pmax(p_obs_po, 1e-15)
ls_po <- -sum(wt_norm * log(p_obs_po), na.rm = TRUE)

# PPOM predicted probs (predict.clm with type="prob" returns P(observed category))
p_obs_ppom <- predict(m_ppom, newdata = trip_count_cz, type = "prob")$fit
p_obs_ppom <- pmax(p_obs_ppom, 1e-15)
ls_ppom <- -sum(wt_norm * log(p_obs_ppom), na.rm = TRUE)

# Calibration for PO model (observed vs predicted weighted bin shares)
obs_bin <- trip_count_cz %>%
  mutate(count_bin = as.character(count_bin)) %>%
  group_by(count_bin) %>%
  summarise(wt = sum(day_weight, na.rm = TRUE), .groups = "drop") %>%
  mutate(p_obs = wt / sum(wt))

p_pred_po <- as.numeric(t(wt_norm) %*% P_po)
calib_po <- data.frame(
  count_bin = bin_levels,
  p_obs  = obs_bin$p_obs[match(bin_levels, obs_bin$count_bin)],
  p_pred = p_pred_po
) %>%
  mutate(abs_diff = abs(p_pred - p_obs))

mae_po <- mean(calib_po$abs_diff)

# Turning point for ordinal model
h_mean <- mean(trip_count_cz$telework_hour, na.rm = TRUE)
h_sd   <- sd(trip_count_cz$telework_hour, na.rm = TRUE)
bz     <- coef(ord_po)["telework_z"]
bz2    <- coef(ord_po)["telework_z2"]
z_star <- -bz / (2 * bz2)
tp_ord <- h_mean + z_star * h_sd

bin_fit <- data.frame(
  Model = c("PO ordered logit (svyolr)", "PPOM (clm, nominal = ~telework_z)"),
  Weighted_log_score = round(c(ls_po, ls_ppom), 5),
  Calibration_MAE = c(round(mae_po, 5), NA),
  N_params = c(length(coef(ord_po)), length(coef(m_ppom))),
  Turning_point = c(round(tp_ord, 2), NA),
  PO_adequate = c(
    ifelse(abs(ls_po - ls_ppom) < 0.001, "Yes", "Marginal"),
    NA
  ),
  stringsAsFactors = FALSE
)

################################################################################
# PANEL C: SENSITIVITY — Linear vs quadratic telework specification
################################################################################

# NB linear telework coefficients
b_lin_tw <- coef(nb_linear)["telework_hour"]
se_lin_tw <- sqrt(diag(cluster_vcov_fn(
  nb_linear,
  trip_count$person_id[as.integer(rownames(model.frame(nb_linear)))]
)))["telework_hour"]

# NB quadratic telework coefficients
cl_q <- trip_count$person_id[as.integer(rownames(model.frame(nb_quad)))]
V_q  <- cluster_vcov_fn(nb_quad, cl_q)
se_q <- sqrt(diag(V_q))

sensitivity_df <- data.frame(
  Check = c(
    "NB: linear telework_hour",
    "NB: quadratic telework_hour",
    "NB: quadratic I(telework_hour^2)",
    "NB: LR test quadratic vs linear (Chi-sq)",
    "NB: AIC improvement (linear - quadratic)",
    "Poisson Pearson overdispersion ratio",
    "NB quadratic: theta (dispersion)",
    "Count model: % zero trips (observed)",
    "Count model: % zero trips (NB predicted)",
    "Ordinal PO: weighted log score",
    "Ordinal PPOM: weighted log score",
    "PO - PPOM log score difference"
  ),
  Estimate = c(
    round(b1_q <- coef(nb_quad)["telework_hour"], 4),
    round(coef(nb_quad)["telework_hour"], 4),
    round(coef(nb_quad)["I(telework_hour^2)"], 4),
    round(as.numeric(-2 * (logLik(nb_linear) - logLik(nb_quad))), 2),
    round(AIC(nb_linear) - AIC(nb_quad), 2),
    round(overdispersion_pear, 2),
    round(nb_quad$theta, 4),
    round(100 * mean(y_all == 0, na.rm = TRUE), 1),
    round(100 * sum(w_all * dnbinom(0, size = theta_q, mu = mu_q), na.rm = TRUE) / tot_w, 1),
    round(ls_po, 5),
    round(ls_ppom, 5),
    round(ls_po - ls_ppom, 5)
  ),
  Robust_SE = c(
    round(se_lin_tw, 4),
    round(se_q["telework_hour"], 4),
    round(se_q["I(telework_hour^2)"], 4),
    NA, NA, NA, NA, NA, NA, NA, NA, NA
  ),
  stringsAsFactors = FALSE
)

################################################################################
# CREATE AND SAVE TABLES
################################################################################

# --- Panel A: Count model comparison ---
tab7a <- gt(count_fit) %>%
  tab_header(
    title = "Panel A: Count model fit comparison",
    subtitle = "Poisson vs NB (linear) vs NB (quadratic telework)"
  ) %>%
  fmt_number(columns = c(AIC, logLik), decimals = 1) %>%
  fmt_number(columns = theta, decimals = 4) %>%
  fmt_number(columns = c(Calibration_MAE), decimals = 5) %>%
  cols_label(
    Model = "Model", AIC = "AIC", logLik = "Log-lik",
    N_params = "k", theta = "Theta",
    Overdispersion = "Overdispersion", Calibration_MAE = "Calib. MAE",
    Turning_point = "h*"
  ) %>%
  sub_missing(missing_text = "--") %>%
  tab_options(table.font.size = px(11))

# --- Panel B: Ordinal model comparison ---
tab7b <- gt(bin_fit) %>%
  tab_header(
    title = "Panel B: Ordinal bin model fit comparison",
    subtitle = "PO (proportional odds) vs PPOM (partial proportional odds)"
  ) %>%
  fmt_number(columns = c(Weighted_log_score, Calibration_MAE), decimals = 5) %>%
  cols_label(
    Model = "Model", Weighted_log_score = "Wt. log score",
    Calibration_MAE = "Calib. MAE", N_params = "k",
    Turning_point = "h*", PO_adequate = "PO adequate?"
  ) %>%
  sub_missing(missing_text = "--") %>%
  tab_options(table.font.size = px(11))

# --- Panel C: Sensitivity checks ---
tab7c <- gt(sensitivity_df) %>%
  tab_header(
    title = "Panel C: Sensitivity checks and diagnostics"
  ) %>%
  cols_label(
    Check = "Diagnostic", Estimate = "Value", Robust_SE = "Robust SE"
  ) %>%
  sub_missing(missing_text = "--") %>%
  tab_source_note(
    "Overdispersion ratio > 1 indicates Poisson inadequacy. Positive AIC improvement favors quadratic. Small |PO - PPOM| log score difference supports proportional odds assumption."
  ) %>%
  tab_options(table.font.size = px(11))

# Save HTML
gt::gtsave(tab7a, "modeling/table-figure/tab7a-count-model-fit.html")
gt::gtsave(tab7b, "modeling/table-figure/tab7b-ordinal-model-fit.html")
gt::gtsave(tab7c, "modeling/table-figure/tab7c-sensitivity-checks.html")


# Print summary
cat("\n=== Panel A: Count model fit ===\n")
print(count_fit)
cat("\n=== Panel B: Ordinal model fit ===\n")
print(bin_fit)
cat("\n=== Panel C: Sensitivity ===\n")
print(sensitivity_df)
