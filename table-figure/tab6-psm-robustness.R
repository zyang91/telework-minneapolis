################################################################################
# Table 6. Propensity-score-matching diagnostics and matched-sample
#           robustness results
################################################################################


library(MASS)
library(sandwich)
library(lmtest)
library(survey)
library(ordinal)
library(gt)

################################################################################
# PANEL A: PSM DIAGNOSTICS - Sample size & balance summary
################################################################################

# Main sample sizes
n_main <- nrow(trip_count)
n_psm  <- nrow(dat_psm)

trip_count <- trip_count %>%
  mutate(telework_binary = ifelse(telework_hour > 0, 1, 0))

dat_psm <- dat_psm %>%
  mutate(telework_binary = as.numeric(as.character(telework_binary)),
         num_people   = as.numeric(num_people),
         num_vehicles = as.numeric(num_vehicles))

# SMD computation
compute_smd_summary <- function(data, treat_var, covariates) {
  smds <- c()
  for (v in covariates) {
    if (is.numeric(data[[v]])) {
      m1 <- mean(data[[v]][data[[treat_var]] == 1], na.rm = TRUE)
      m0 <- mean(data[[v]][data[[treat_var]] == 0], na.rm = TRUE)
      s1 <- sd(data[[v]][data[[treat_var]] == 1], na.rm = TRUE)
      s0 <- sd(data[[v]][data[[treat_var]] == 0], na.rm = TRUE)
      pooled_sd <- sqrt((s1^2 + s0^2) / 2)
      smds <- c(smds, abs((m1 - m0) / pooled_sd))
    } else {
      levs <- unique(data[[v]])
      for (lv in levs) {
        m1 <- mean(data[[v]][data[[treat_var]] == 1] == lv, na.rm = TRUE)
        m0 <- mean(data[[v]][data[[treat_var]] == 0] == lv, na.rm = TRUE)
        pooled_sd <- sqrt((m1 * (1 - m1) + m0 * (1 - m0)) / 2)
        if (pooled_sd > 0) smds <- c(smds, abs((m1 - m0) / pooled_sd))
      }
    }
  }
  smds
}

balance_vars <- c("age", "gender", "education", "num_people", "num_vehicles",
                  "income_broad", "thrive_community_type", "travel_dow",
                  "travel_date_season")

smd_before <- compute_smd_summary(trip_count, "telework_binary", balance_vars)
smd_after  <- compute_smd_summary(dat_psm, "telework_binary", balance_vars)

diag_df <- data.frame(
  Metric = c("N (person-days)", "N (unique persons)",
             "Teleworkers (N)", "Non-teleworkers (N)",
             "Mean |SMD| before matching", "Max |SMD| before matching",
             "Mean |SMD| after matching", "Max |SMD| after matching",
             "Covariates with |SMD| > 0.1 before", "Covariates with |SMD| > 0.1 after"),
  Value = c(
    format(n_psm, big.mark = ","),
    format(length(unique(dat_psm$person_id)), big.mark = ","),
    format(sum(dat_psm$telework_binary == 1, na.rm = TRUE), big.mark = ","),
    format(sum(dat_psm$telework_binary == 0, na.rm = TRUE), big.mark = ","),
    sprintf("%.4f", mean(smd_before, na.rm = TRUE)),
    sprintf("%.4f", max(smd_before, na.rm = TRUE)),
    sprintf("%.4f", mean(smd_after, na.rm = TRUE)),
    sprintf("%.4f", max(smd_after, na.rm = TRUE)),
    as.character(sum(smd_before > 0.1)),
    as.character(sum(smd_after > 0.1))
  )
)

################################################################################
# PANEL B: PSM ROBUSTNESS - Telework coefficients across models
################################################################################

# --- B1: PSM Ordered logit (bins) ---
dat_psm <- dat_psm %>%
  mutate(
    telework_c  = telework_hour - mean(telework_hour, na.rm = TRUE),
    telework_z  = as.numeric(scale(telework_c)),
    telework_z2 = telework_z^2
  )

ord_psm <- clm(
  count_bin ~ travel_dow + age + gender + education +
    telework_z + telework_z2 +
    num_people + num_vehicles +
    income_broad + thrive_community_type + travel_date_season + year,
  data = dat_psm,
  link = "logit",
  Hess = TRUE
)

ord_cf <- summary(ord_psm)$coefficients
ord_tele <- ord_cf[grepl("telework", rownames(ord_cf)), , drop = FALSE]

# --- B2: PSM NB (total trips) ---
nb_psm <- glm.nb(
  count ~ travel_dow + age + gender + education +
    telework_hour + I(telework_hour^2) +
    num_people + num_vehicles +
    income_broad + thrive_community_type + travel_date_season + year,
  data = dat_psm
)

mf_nb_psm <- model.frame(nb_psm)
idx_nb_psm <- as.integer(rownames(mf_nb_psm))
cl_nb_psm <- dat_psm$person_id[idx_nb_psm]
V_nb_psm <- cluster_vcov_fn(nb_psm, cl_nb_psm)

b_nb  <- coef(nb_psm)
se_nb <- sqrt(diag(V_nb_psm))
nb_tele_idx <- grepl("telework", names(b_nb))

# --- B3: PSM Work logistic ---
dat_psm <- dat_psm %>%
  mutate(work_day_bin = ifelse(work_day == "yes", 1, 0))

m_work_psm <- glm(
  work_day_bin ~ travel_dow + age + gender + education +
    telework_hour + I(telework_hour^2) +
    num_people + num_vehicles +
    income_broad + thrive_community_type + travel_date_season + year,
  data = dat_psm,
  family = binomial()
)

work_cf <- summary(m_work_psm)$coefficients
work_tele <- work_cf[grepl("telework", rownames(work_cf)), , drop = FALSE]

# --- B4: PSM Non-work NB ---
dat_psm$non_work_trip[is.na(dat_psm$non_work_trip)] <- 0

nb_nw_psm <- glm.nb(
  non_work_trip ~ travel_dow + age + gender + education +
    telework_hour + I(telework_hour^2) +
    num_people + num_vehicles +
    income_broad + thrive_community_type + travel_date_season + year,
  data = dat_psm
)

mf_nw_psm <- model.frame(nb_nw_psm)
idx_nw_psm <- as.integer(rownames(mf_nw_psm))
cl_nw_psm <- dat_psm$person_id[idx_nw_psm]
V_nw_psm <- cluster_vcov_fn(nb_nw_psm, cl_nw_psm)

b_nw_psm  <- coef(nb_nw_psm)
se_nw_psm <- sqrt(diag(V_nw_psm))

# --- Combine telework results ---
robustness_df <- data.frame(
  Model = c(
    "Ordered logit (bins) - telework_z",
    "Ordered logit (bins) - telework_z2",
    "NB total trips - telework_hour",
    "NB total trips - I(telework_hour^2)",
    "Work logistic - telework_hour",
    "Work logistic - I(telework_hour^2)",
    "NB non-work - telework_hour",
    "NB non-work - I(telework_hour^2)"
  ),
  Estimate = round(c(
    ord_tele[, "Estimate"],
    b_nb[nb_tele_idx],
    work_tele[, "Estimate"],
    b_nw_psm[grepl("telework", names(b_nw_psm))]
  ), 4),
  SE = round(c(
    ord_tele[, "Std. Error"],
    se_nb[nb_tele_idx],
    work_tele[, "Std. Error"],
    se_nw_psm[grepl("telework", names(se_nw_psm))]
  ), 4),
  p_value = c(
    ord_tele[, "Pr(>|z|)"],
    2 * pnorm(-abs(b_nb[nb_tele_idx] / se_nb[nb_tele_idx])),
    work_tele[, "Pr(>|z|)"],
    2 * pnorm(-abs(b_nw_psm[grepl("telework", names(b_nw_psm))] /
                     se_nw_psm[grepl("telework", names(se_nw_psm))]))
  ),
  stringsAsFactors = FALSE,
  row.names = NULL
)

robustness_df$sig <- case_when(
  robustness_df$p_value < 0.001 ~ "***",
  robustness_df$p_value < 0.01  ~ "**",
  robustness_df$p_value < 0.05  ~ "*",
  robustness_df$p_value < 0.1   ~ ".",
  TRUE                           ~ ""
)
robustness_df$p_value <- ifelse(robustness_df$p_value < 0.001, "<0.001",
                                sprintf("%.3f", robustness_df$p_value))

################################################################################
# SAVE TABLES
################################################################################

tab6a <- gt(diag_df) %>%
  tab_header(title = "Panel A: PSM diagnostics") %>%
  cols_label(Metric = "Diagnostic", Value = "Value") %>%
  tab_options(table.font.size = px(11))

tab6b <- gt(robustness_df) %>%
  tab_header(
    title = "Panel B: Telework coefficients in PSM matched sample",
    subtitle = "Robustness comparison across all model specifications"
  ) %>%
  fmt_number(columns = c(Estimate, SE), decimals = 4) %>%
  cols_label(
    Model = "Model / Variable", Estimate = "Coef.", SE = "SE",
    p_value = "p", sig = ""
  ) %>%
  tab_source_note("Signif. codes: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1") %>%
  tab_options(table.font.size = px(11))

gt::gtsave(tab6a, "modeling/table-figure/tab6a-psm-diagnostics.html")
gt::gtsave(tab6b, "modeling/table-figure/tab6b-psm-robustness.html")
