################################################################################
# Table 4. Negative binomial estimates for total daily trip counts
################################################################################


library(MASS)
library(sandwich)
library(lmtest)
library(gt)

# --- 1. Fit quadratic NB model ---
nb_model_q <- glm.nb(
  count ~ travel_dow + age + gender + education +
    telework_hour + I(telework_hour^2) +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  data = trip_count,
  weights = trip_count$day_weight
)

# --- 2. Cluster-robust SEs ---
mf_q  <- model.frame(nb_model_q)
idx_q <- as.integer(rownames(mf_q))
cluster_q <- trip_count$person_id[idx_q]

V_cl_q <- cluster_vcov_fn(nb_model_q, cluster_q)

b_q  <- coef(nb_model_q)
se_q <- sqrt(diag(V_cl_q))
z_q  <- b_q / se_q
p_q  <- 2 * pnorm(abs(z_q), lower.tail = FALSE)

irr_tab <- data.frame(
  Variable  = names(b_q),
  Estimate  = round(as.numeric(b_q), 4),
  SE_robust = round(as.numeric(se_q), 4),
  z         = round(as.numeric(z_q), 3),
  p_value   = as.numeric(p_q),
  IRR       = round(exp(as.numeric(b_q)), 4),
  IRR_lo    = round(exp(as.numeric(b_q) - 1.96 * as.numeric(se_q)), 4),
  IRR_hi    = round(exp(as.numeric(b_q) + 1.96 * as.numeric(se_q)), 4),
  stringsAsFactors = FALSE,
  row.names = NULL
)

# Add significance stars
irr_tab$sig <- case_when(
  irr_tab$p_value < 0.001 ~ "***",
  irr_tab$p_value < 0.01  ~ "**",
  irr_tab$p_value < 0.05  ~ "*",
  irr_tab$p_value < 0.1   ~ ".",
  TRUE                     ~ ""
)

# Format p-values
irr_tab$p_value <- ifelse(irr_tab$p_value < 0.001, "<0.001",
                          sprintf("%.3f", irr_tab$p_value))

# --- 3. Turning point ---
b1 <- coef(nb_model_q)["telework_hour"]
b2 <- coef(nb_model_q)["I(telework_hour^2)"]
tp <- -b1 / (2 * b2)

# --- 4. Create gt table ---
tab4 <- gt(irr_tab) %>%
  tab_header(
    title = "Negative binomial estimates for total daily trip counts",
    subtitle = "Quadratic telework specification; cluster-robust SEs (person_id)"
  ) %>%
  fmt_number(columns = c(Estimate, SE_robust, IRR, IRR_lo, IRR_hi), decimals = 4) %>%
  fmt_number(columns = z, decimals = 3) %>%
  cols_label(
    Variable  = "Variable",
    Estimate  = "Coef.",
    SE_robust = "Robust SE",
    z         = "z",
    p_value   = "p",
    IRR       = "IRR",
    IRR_lo    = "IRR 2.5%",
    IRR_hi    = "IRR 97.5%",
    sig       = ""
  ) %>%
  cols_width(Variable ~ px(250)) %>%
  tab_source_note(
    sprintf("N = %d observations, %d person-clusters. Theta = %.4f. Turning point h* = %.2f hrs.",
            nobs(nb_model_q), nlevels(as.factor(cluster_q)),
            nb_model_q$theta, tp)
  ) %>%
  tab_source_note("Signif. codes: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1") %>%
  tab_options(table.font.size = px(11))

# --- 5. Save ---
gt::gtsave(tab4, "modeling/table-figure/tab4-nb-estimates.html")

