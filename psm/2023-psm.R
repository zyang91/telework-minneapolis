################################################################################
# 2023 Propensity Score Matching (PSM) Script
################################################################################
#
# Purpose:
#   This script implements propensity score matching for the 2023 Travel Behavior
#   Inventory (TBI) survey data. The goal is to create a balanced comparison
#   between person-days with and without working from home (WFH), so that
#   differences in travel behavior can be attributed to telework rather than
#   pre-existing differences between teleworkers and non-teleworkers.
#
# Treatment Variable:
#   - telework_binary: 1 = worked from home (>6 hours on that day), 0 = did not
#
# Analysis Unit:
#   - Person-day (each row is one person on one survey day)
#
# Matching Covariates (from the main model in modeling/model1/trip_count.R):
#   - Demographics: age, gender, education
#   - Employment: telework_freq, job_type, hours_work, work_mode
#   - Household: num_people, num_vehicles, res_type, income_broad
#   - Location: thrive_community_type
#   - Temporal: travel_dow, travel_date_season
#
# Steps:
#   1. Load and prepare the 2023 person-day data
#   2. Harmonize categorical variables (same rules as main model)
#   3. Check treatment distribution
#   4. Estimate propensity scores using logistic regression
#   5. Perform nearest-neighbor 1:1 matching without replacement
#   6. Evaluate covariate balance (before vs. after matching)
#   7. Generate love plot to visualize balance
#   8. Report evaluation statistics
#
# Input:
#   - YOUR TBI FILE PATH HERE (e.g., "data/2023-clean/trip_count_clean_2023.csv")
#
# Output:
#   - Matched dataset stored in `matched_data`
#   - Balance table printed to console
#   - Love plot saved as output
#
# Author: Zhanchao Yang, University of Pennsylvania
# Last Modified: June 2025
# R Version: 4.0+
#
################################################################################


# ── Libraries ─────────────────────────────────────────────────────────────────

library(tidyverse)   # data manipulation and ggplot2
library(MatchIt)     # propensity score matching
library(cobalt)      # love plot and balance statistics
options(scipen = 999)


# ── Step 1: Load Data ─────────────────────────────────────────────────────────
# The cleaned 2023 person-day dataset already has all the variables we need.

data_2023 <- read.csv("YOUR FILE PATH")



# ── Step 2: Harmonize Categorical Variables ───────────────────────────────────
# Apply the same recoding rules used in the main model (modeling/model1/trip_count.R)
# so that the PSM uses consistent variable definitions.

# Parse household size and vehicle count from strings to numbers
data_2023 <- data_2023 %>%
  mutate(
    num_people   = parse_number(as.character(num_people)),
    num_vehicles = parse_number(as.character(num_vehicles))
  )

## other categorical variables haramonize (education, job_type, res_type, income_broad, etc.)

# ── Step 3: Check Treatment Distribution ──────────────────────────────────────
# The treatment is whether the person worked from home on that day.
# telework_binary = 1  →  WFH day (>6 hours working from home)
# telework_binary = 0  →  non-WFH day

print(table(data_2023$telework_binary))

cat("\nTreatment proportion:", round(mean(data_2023$telework_binary), 3), "\n")


# ── Step 4: Prepare Data for Matching ─────────────────────────────────────────
# Select only the variables needed for matching.
# We exclude outcome variables (count, non_work_trip, work_day, etc.) so
# that the propensity score only reflects pre-treatment characteristics.

psm_data <- data_2023 %>%
  dplyr::select(
    # Treatment
    telework_binary,
    # Person-level demographics
    age, gender, education,
    # Employment characteristics (predict whether someone WFH on a given day)
    telework_freq, job_type, hours_work, work_mode,
    # Household characteristics
    num_people, num_vehicles, res_type, income_broad,
    # Location
    thrive_community_type,
    # Temporal context
    travel_dow, travel_date_season,
    # IDs for reference (not used in model)
    person_id, day_id, day_weight
  ) %>%
  # Drop rows with any missing values in the matching covariates
  # (explicitly list each variable so the behavior is clear)
  drop_na(age, gender, education,
          telework_freq, job_type, hours_work, work_mode,
          num_people, num_vehicles, res_type, income_broad,
          thrive_community_type, travel_dow, travel_date_season)

cat("\nPerson-days available for matching:", nrow(psm_data), "\n")
cat("  WFH days (treated):     ", sum(psm_data$telework_binary == 1), "\n")
cat("  Non-WFH days (control): ", sum(psm_data$telework_binary == 0), "\n")


# ── Step 5: Propensity Score Matching Using MatchIt ───────────────────────────
# We use nearest-neighbor 1:1 matching (each WFH day is matched to the most
# similar non-WFH day based on estimated propensity scores).
#
# The propensity score is the predicted probability of WFH on a given day,
# estimated by logistic regression using all the covariates above.
#
# method = "nearest"   → nearest-neighbor matching
# ratio  = 1           → 1:1 matching (one control per treated unit)
# replace = FALSE      → each control unit can only be used once

set.seed(123)  # for reproducibility

match_out <- matchit(
  telework_binary ~
    age + gender + education +
    hours_work +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_dow + travel_date_season,
  data    = psm_data,
  method  = "nearest",
  distance = "logit",   # estimate propensity score via logistic regression
  ratio   = 1,
  replace = FALSE,
  caliper  = 0.2,
  std.caliper = TRUE
)

# Print a brief summary of the matching procedure
print(summary(match_out, standardize = TRUE))


# ── Step 6: Extract the Matched Dataset ───────────────────────────────────────
# match.data() returns only the matched observations with an added column
# `weights` that reflects matching weights (1 for matched, 0 for discarded).

matched_data <- match.data(match_out)

cat("\nMatched dataset size:", nrow(matched_data), "person-days\n")
cat("  WFH days (treated):     ", sum(matched_data$telework_binary == 1), "\n")
cat("  Non-WFH days (control): ", sum(matched_data$telework_binary == 0), "\n")


# ── Step 7: Balance Evaluation ────────────────────────────────────────────────
# We compare the standardized mean differences (SMD) before and after matching.
# A well-matched dataset should have SMD close to 0 for all covariates.
# A common rule of thumb: SMD < 0.1 indicates good balance.

bal_table <- bal.tab(
  match_out,
  stats    = c("mean.diffs", "variance.ratios"),
  thresholds = c(m = 0.1, v = 2),  # SMD < 0.1 and variance ratio < 2
  un       = TRUE   # also show balance BEFORE matching
)
print(bal_table)


# ── Step 8: Love Plot ─────────────────────────────────────────────────────────
# A love plot visualizes the standardized mean difference for each covariate
# before and after matching. Points near the vertical line at 0 indicate
# good balance.

love.plot(
  match_out,
  stats        = "mean.diffs",
  threshold    = 0.1,           # dashed line showing the "good balance" threshold
  binary       = "std",         # standardize binary variables too
  abs          = TRUE,          # show absolute SMD
  var.order    = "unadjusted",  # order variables by pre-match SMD
  shapes       = c("circle filled", "triangle filled"),
  colors       = c("#E87722", "#003865"),
  title        = "Covariate Balance Before and After Matching (2023 Data)",
  sample.names = c("Before Matching", "After Matching")
)


# ── Step 9: Additional Evaluation Statistics ──────────────────────────────────
# Report key balance statistics as a simple data frame.

bal_stats <- bal.tab(
  match_out,
  stats = c("mean.diffs", "variance.ratios"),
  un    = TRUE
)

# Extract the balance table as a data frame
balance_df <- bal_stats$Balance

cat("\n── Variables with |SMD| > 0.1 after matching (needs review) ─────────\n")
# Standardized mean difference after matching is in column "Diff.Adj"
if ("Diff.Adj" %in% names(balance_df)) {
  poor_balance <- balance_df %>%
    rownames_to_column("variable") %>%
    filter(abs(Diff.Adj) > 0.1)
  if (nrow(poor_balance) == 0) {
    cat("All covariates have |SMD| <= 0.1 after matching. Balance is good.\n")
  } else {
    print(poor_balance[, c("variable", "Diff.Un", "Diff.Adj")])
  }
}

# Summary of propensity score distribution by treatment group
psm_data$propensity_score <- match_out$distance

cat("\n── Propensity Score Distribution ─────────────────────────────────────\n")
ps_summary <- psm_data %>%
  group_by(telework_binary) %>%
  summarise(
    n        = n(),
    mean_ps  = round(mean(propensity_score), 4),
    sd_ps    = round(sd(propensity_score), 4),
    min_ps   = round(min(propensity_score), 4),
    max_ps   = round(max(propensity_score), 4),
    .groups  = "drop"
  ) %>%
  mutate(group = ifelse(telework_binary == 1, "WFH day", "Non-WFH day"))

print(ps_summary)

# Propensity score overlap plot (check common support)
ggplot(psm_data, aes(x = propensity_score,
                     fill  = factor(telework_binary),
                     color = factor(telework_binary))) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(
    values = c("0" = "#003865", "1" = "#E87722"),
    labels = c("0" = "Non-WFH day", "1" = "WFH day")
  ) +
  scale_color_manual(
    values = c("0" = "#003865", "1" = "#E87722"),
    labels = c("0" = "Non-WFH day", "1" = "WFH day")
  ) +
  labs(
    title   = "Propensity Score Distribution by Treatment Group (2023)",
    x       = "Propensity Score (Estimated Probability of WFH)",
    y       = "Density",
    fill    = "Group",
    color   = "Group",
    caption = "Good overlap between the two distributions indicates common support."
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")


# export matched_data
# write.csv(matched_data, "data/2023-clean/psm_matched_data_2023.csv", row.names = FALSE)
