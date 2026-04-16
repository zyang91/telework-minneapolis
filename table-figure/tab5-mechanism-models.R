################################################################################
# Table 5. Mechanism model estimates for work-related travel and non-work trips
################################################################################

library(MASS)
library(sandwich)
library(lmtest)
library(survey)
library(gt)

################################################################################
# PANEL A: Work-travel logistic regression (substitution channel)
################################################################################

trip_count <- trip_count %>%
  mutate(
    work_day_bin = ifelse(work_day == "yes", 1, 0),
    telework_c   = telework_hour - mean(telework_hour, na.rm = TRUE),
    telework_c2  = telework_c^2
  )

des <- svydesign(ids = ~person_id, weights = ~day_weight, data = trip_count)

m_work_c <- svyglm(
  work_day_bin ~ travel_dow + age + gender + education +
    telework_c + telework_c2 +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  design = des,
  family = quasibinomial()
)

# Extract coefficients
sm_work <- summary(m_work_c)
ct_work <- sm_work$coefficients

work_df <- data.frame(
  Variable = rownames(ct_work),
  Estimate = round(ct_work[, "Estimate"], 4),
  SE       = round(ct_work[, "Std. Error"], 4),
  z_value  = round(ct_work[, "t value"], 3),
  p_value  = ct_work[, "Pr(>|t|)"],
  OR       = round(exp(ct_work[, "Estimate"]), 4),
  OR_lo    = round(exp(ct_work[, "Estimate"] - 1.96 * ct_work[, "Std. Error"]), 4),
  OR_hi    = round(exp(ct_work[, "Estimate"] + 1.96 * ct_work[, "Std. Error"]), 4),
  stringsAsFactors = FALSE,
  row.names = NULL
)
work_df$sig <- case_when(
  work_df$p_value < 0.001 ~ "***",
  work_df$p_value < 0.01  ~ "**",
  work_df$p_value < 0.05  ~ "*",
  work_df$p_value < 0.1   ~ ".",
  TRUE                     ~ ""
)
work_df$p_value <- ifelse(work_df$p_value < 0.001, "<0.001",
                          sprintf("%.3f", work_df$p_value))
work_df$model <- "Work-travel logistic"

################################################################################
# PANEL B: Non-work trips NB (supplementation channel)
################################################################################

trip_count$non_work_trip[is.na(trip_count$non_work_trip)] <- 0

nb_nonwork <- glm.nb(
  non_work_trip ~ travel_dow + age + gender + education +
    telework_hour + I(telework_hour^2) +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  data = trip_count,
  weights = trip_count$day_weight
)

# Cluster-robust SEs
mf_nw  <- model.frame(nb_nonwork)
idx_nw <- as.integer(rownames(mf_nw))
cluster_nw <- trip_count$person_id[idx_nw]
V_cl_nw <- cluster_vcov_fn(nb_nonwork, cluster_nw)

b_nw  <- coef(nb_nonwork)
se_nw <- sqrt(diag(V_cl_nw))
z_nw  <- b_nw / se_nw
p_nw  <- 2 * pnorm(abs(z_nw), lower.tail = FALSE)

nw_df <- data.frame(
  Variable  = names(b_nw),
  Estimate  = round(as.numeric(b_nw), 4),
  SE        = round(as.numeric(se_nw), 4),
  z_value   = round(as.numeric(z_nw), 3),
  p_value   = as.numeric(p_nw),
  IRR       = round(exp(as.numeric(b_nw)), 4),
  IRR_lo    = round(exp(as.numeric(b_nw) - 1.96 * as.numeric(se_nw)), 4),
  IRR_hi    = round(exp(as.numeric(b_nw) + 1.96 * as.numeric(se_nw)), 4),
  stringsAsFactors = FALSE,
  row.names = NULL
)
nw_df$sig <- case_when(
  nw_df$p_value < 0.001 ~ "***",
  nw_df$p_value < 0.01  ~ "**",
  nw_df$p_value < 0.05  ~ "*",
  nw_df$p_value < 0.1   ~ ".",
  TRUE                   ~ ""
)
nw_df$p_value <- ifelse(nw_df$p_value < 0.001, "<0.001",
                        sprintf("%.3f", nw_df$p_value))
nw_df$model <- "Non-work NB"

################################################################################
# COMBINED TABLE (side-by-side or stacked)
################################################################################

# Save separate gt tables

tab5a <- gt(work_df %>% select(-model)) %>%
  tab_header(
    title = "Panel A: Work-travel logistic regression (substitution channel)",
    subtitle = "Survey-weighted quasibinomial (svyglm); person-cluster design; centered telework"
  ) %>%
  fmt_number(columns = c(Estimate, SE, OR, OR_lo, OR_hi), decimals = 4) %>%
  fmt_number(columns = z_value, decimals = 3) %>%
  cols_label(
    Variable = "Variable", Estimate = "Coef.", SE = "SE", z_value = "t",
    p_value = "p", OR = "OR", OR_lo = "OR 2.5%", OR_hi = "OR 97.5%", sig = ""
  ) %>%
  tab_source_note("Signif. codes: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1") %>%
  tab_options(table.font.size = px(11))

tab5b <- gt(nw_df %>% select(-model)) %>%
  tab_header(
    title = "Panel B: Non-work trips negative binomial (supplementation channel)",
    subtitle = "Quadratic telework specification; cluster-robust SEs (person_id)"
  ) %>%
  fmt_number(columns = c(Estimate, SE, IRR, IRR_lo, IRR_hi), decimals = 4) %>%
  fmt_number(columns = z_value, decimals = 3) %>%
  cols_label(
    Variable = "Variable", Estimate = "Coef.", SE = "Robust SE", z_value = "z",
    p_value = "p", IRR = "IRR", IRR_lo = "IRR 2.5%", IRR_hi = "IRR 97.5%", sig = ""
  ) %>%
  tab_source_note(sprintf("N = %d, theta = %.4f", nobs(nb_nonwork), nb_nonwork$theta)) %>%
  tab_source_note("Signif. codes: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1") %>%
  tab_options(table.font.size = px(11))

# Save
gt::gtsave(tab5a, "modeling/table-figure/tab5a-work-logistic.html")
gt::gtsave(tab5b, "modeling/table-figure/tab5b-nonwork-nb.html")

# Turning points
tw_mean <- mean(trip_count$telework_hour, na.rm = TRUE)
b_c1 <- coef(m_work_c)["telework_c"]
b_c2 <- coef(m_work_c)["telework_c2"]
h_star_work <- -b_c1 / (2 * b_c2) + tw_mean

b_nw1 <- coef(nb_nonwork)["telework_hour"]
b_nw2 <- coef(nb_nonwork)["I(telework_hour^2)"]
h_star_nw <- -b_nw1 / (2 * b_nw2)

cat(sprintf("\nWork-travel turning point h*: %.2f hrs\n", h_star_work))
cat(sprintf("Non-work turning point h*:   %.2f hrs\n", h_star_nw))

