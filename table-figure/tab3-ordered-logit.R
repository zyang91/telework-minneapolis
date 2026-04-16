################################################################################
# Table 3. Ordered logit estimates for daily trip-count bins
################################################################################


library(survey)
library(gt)

# --- 1. Prepare data with standardized telework ---
trip_count_cz <- trip_count %>%
  mutate(
    telework_c  = telework_hour - mean(telework_hour, na.rm = TRUE),
    telework_z  = as.numeric(scale(telework_c)),
    telework_z2 = telework_z^2
  )

# --- 2. Fit survey-weighted ordered logit ---
svy_design <- svydesign(ids = ~person_id, weights = ~day_weight, data = trip_count_cz)

ord_mod_z <- svyolr(
  count_bin ~ travel_dow + age + gender + education +
    telework_z + telework_z2 +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  design = svy_design
)

# --- 3. Extract coefficients ---
smry_coef <- summary(ord_mod_z)$coefficients

b  <- smry_coef[, "Value"]
se <- smry_coef[, "Std. Error"]
z  <- smry_coef[, "t value"]
p  <- 2 * pnorm(-abs(z))

coef_df <- data.frame(
  Variable = rownames(smry_coef),
  Estimate = round(as.numeric(b), 4),
  SE       = round(as.numeric(se), 4),
  z_value  = round(as.numeric(z), 3),
  p_value  = as.numeric(p),
  OR       = round(exp(as.numeric(b)), 4),
  OR_lo    = round(exp(as.numeric(b) - 1.96 * as.numeric(se)), 4),
  OR_hi    = round(exp(as.numeric(b) + 1.96 * as.numeric(se)), 4),
  stringsAsFactors = FALSE,
  row.names = NULL
)

# Separate thresholds and slopes
is_thresh <- grepl("\\|", coef_df$Variable)
thresh_df <- coef_df[is_thresh, ]
slope_df  <- coef_df[!is_thresh, ]

# Add significance stars
add_stars <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*",
            p < 0.1 ~ ".", TRUE ~ "")
}
slope_df$sig <- add_stars(slope_df$p_value)

# Format p-values
slope_df$p_value <- ifelse(slope_df$p_value < 0.001, "<0.001",
                           sprintf("%.3f", slope_df$p_value))

# --- 4. Create gt table ---
tab3 <- gt(slope_df) %>%
  tab_header(
    title = "Ordered logit estimates for daily trip-count bins",
    subtitle = "Survey-weighted (svyolr); person-cluster design; standardized telework"
  ) %>%
  fmt_number(columns = c(Estimate, SE, OR, OR_lo, OR_hi), decimals = 4) %>%
  fmt_number(columns = z_value, decimals = 3) %>%
  cols_label(
    Variable = "Variable",
    Estimate = "Coef.",
    SE       = "SE",
    z_value  = "z",
    p_value  = "p",
    OR       = "OR",
    OR_lo    = "OR 2.5%",
    OR_hi    = "OR 97.5%",
    sig      = ""
  ) %>%
  cols_width(Variable ~ px(250)) %>%
  tab_source_note(
    "Bins: 0 | 1-2 | 3-4 | 5-6 | 7+. Positive coefficient = higher odds of being in a higher trip bin."
  ) %>%
  tab_source_note(
    "Signif. codes: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1"
  ) %>%
  tab_options(table.font.size = px(11))

# Threshold table
tab3_thresh <- gt(thresh_df %>% select(Variable, Estimate, SE, z_value)) %>%
  tab_header(title = "Threshold parameters") %>%
  fmt_number(columns = c(Estimate, SE), decimals = 4) %>%
  fmt_number(columns = z_value, decimals = 3)

# --- 5. Model summary stats ---
h_mean <- mean(trip_count_cz$telework_hour, na.rm = TRUE)
h_sd   <- sd(trip_count_cz$telework_hour, na.rm = TRUE)
b_z    <- coef(ord_mod_z)["telework_z"]
b_z2   <- coef(ord_mod_z)["telework_z2"]
z_star <- -b_z / (2 * b_z2)
h_star <- h_mean + z_star * h_sd

cat(sprintf("N observations: %d\n", nrow(trip_count_cz)))
cat(sprintf("N clusters (persons): %d\n", length(unique(trip_count_cz$person_id))))
cat(sprintf("Telework scale: z = (h - %.3f) / %.3f\n", h_mean, h_sd))
cat(sprintf("Turning point h*: %.2f hours\n", h_star))

# --- 6. Save ---
gt::gtsave(tab3, "modeling/table-figure/tab3-ordered-logit.html")

