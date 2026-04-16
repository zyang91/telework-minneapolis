################################################################################
# Table 8. Consolidated turning-point comparison across all model specifications
#
# Shows h* (telework hours at which the nonlinear effect peaks/troughs) for
# every model in the paper, plus the PSM robustness variants.
################################################################################

library(MASS)
library(sandwich)
library(survey)
library(ordinal)
library(gt)

################################################################################
# 1. MAIN-SAMPLE MODELS
################################################################################

# --- A. NB quadratic (total trips) ---
nb_quad <- glm.nb(
  count ~ travel_dow + age + gender + education +
    telework_hour + I(telework_hour^2) +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  data = trip_count,
  weights = trip_count$day_weight
)
tp_nb <- -coef(nb_quad)["telework_hour"] / (2 * coef(nb_quad)["I(telework_hour^2)"])

# --- B. Ordered logit (trip bins, svyolr) ---
trip_count_cz <- trip_count %>%
  mutate(
    telework_c  = telework_hour - mean(telework_hour, na.rm = TRUE),
    telework_z  = as.numeric(scale(telework_c)),
    telework_z2 = telework_z^2
  )
svy_des <- svydesign(ids = ~person_id, weights = ~day_weight, data = trip_count_cz)
ord_po <- svyolr(
  count_bin ~ travel_dow + age + gender + education +
    telework_z + telework_z2 +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  design = svy_des
)
h_mean <- mean(trip_count_cz$telework_hour, na.rm = TRUE)
h_sd   <- sd(trip_count_cz$telework_hour, na.rm = TRUE)
z_star_po <- -coef(ord_po)["telework_z"] / (2 * coef(ord_po)["telework_z2"])
tp_ord <- h_mean + z_star_po * h_sd

# --- C. Non-work NB (supplementation) ---
trip_count$non_work_trip[is.na(trip_count$non_work_trip)] <- 0
nb_nw <- glm.nb(
  non_work_trip ~ travel_dow + age + gender + education +
    telework_hour + I(telework_hour^2) +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  data = trip_count,
  weights = trip_count$day_weight
)
tp_nw <- -coef(nb_nw)["telework_hour"] / (2 * coef(nb_nw)["I(telework_hour^2)"])

################################################################################
# 2. PSM ROBUSTNESS MODELS
################################################################################

# --- D. PSM NB quadratic (total trips) ---
nb_psm <- glm.nb(
  count ~ travel_dow + age + gender + education +
    telework_hour + I(telework_hour^2) +
    num_people + num_vehicles +
    income_broad + thrive_community_type + travel_date_season + year,
  data = dat_psm
)
tp_psm_nb <- -coef(nb_psm)["telework_hour"] / (2 * coef(nb_psm)["I(telework_hour^2)"])

# --- E. PSM ordered logit (trip bins) ---
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
h_mean_psm <- mean(dat_psm$telework_hour, na.rm = TRUE)
h_sd_psm   <- sd(dat_psm$telework_hour, na.rm = TRUE)
z_star_psm <- -coef(ord_psm)["telework_z"] / (2 * coef(ord_psm)["telework_z2"])
tp_psm_ord <- h_mean_psm + z_star_psm * h_sd_psm

# --- F. PSM non-work NB ---
dat_psm$non_work_trip[is.na(dat_psm$non_work_trip)] <- 0
nb_nw_psm <- glm.nb(
  non_work_trip ~ travel_dow + age + gender + education +
    telework_hour + I(telework_hour^2) +
    num_people + num_vehicles +
    income_broad + thrive_community_type + travel_date_season + year,
  data = dat_psm
)
tp_psm_nw <- -coef(nb_nw_psm)["telework_hour"] /
  (2 * coef(nb_nw_psm)["I(telework_hour^2)"])

################################################################################
# 3. BUILD TABLE
################################################################################

tp_table <- data.frame(
  Model = c(
    "NB total trips (quadratic)",
    "Ordered logit trip bins (PO)",
    "Non-work trips NB",
    "PSM: NB total trips",
    "PSM: Ordered logit bins",
    "PSM: Non-work trips NB"
  ),
  Sample = c(rep("Full sample", 3), rep("PSM matched", 3)),
  Outcome = c("Total trips", "Trip bins", "Non-work trips",
              "Total trips", "Trip bins", "Non-work trips"),
  Turning_point_hrs = round(c(tp_nb, tp_ord, tp_nw,
                              tp_psm_nb, tp_psm_ord, tp_psm_nw), 2),
  Shape = ifelse(
    c(coef(nb_quad)["I(telework_hour^2)"],
      coef(ord_po)["telework_z2"],
      coef(nb_nw)["I(telework_hour^2)"],
      coef(nb_psm)["I(telework_hour^2)"],
      coef(ord_psm)["telework_z2"],
      coef(nb_nw_psm)["I(telework_hour^2)"]) < 0,
    "Inverted-U", "U-shape"
  ),
  stringsAsFactors = FALSE
)

# Summary stats
cat(sprintf("Turning point range: [%.2f, %.2f] hours\n",
    min(tp_table$Turning_point_hrs), max(tp_table$Turning_point_hrs)))
cat(sprintf("Mean turning point:  %.2f hours\n",
    mean(tp_table$Turning_point_hrs)))
cat(sprintf("SD of turning points: %.2f hours\n",
    sd(tp_table$Turning_point_hrs)))

################################################################################
# 4. SAVE
################################################################################

tab8 <- gt(tp_table) %>%
  tab_header(
    title = "Turning-point estimates across all model specifications",
    subtitle = "h* = telework hours at which the quadratic effect peaks or troughs"
  ) %>%
  cols_label(
    Model = "Model",
    Sample = "Sample",
    Outcome = "Outcome",
    Turning_point_hrs = "h* (hours)",
    Shape = "Curve shape"
  ) %>%
  tab_row_group(label = "Full sample", rows = 1:4) %>%
  tab_row_group(label = "PSM matched sample", rows = 5:8) %>%
  tab_source_note(
    sprintf("Turning point range: [%.2f, %.2f] hrs. Mean = %.2f hrs.",
            min(tp_table$Turning_point_hrs),
            max(tp_table$Turning_point_hrs),
            mean(tp_table$Turning_point_hrs))
  ) %>%
  tab_options(table.font.size = px(11))

gt::gtsave(tab8, "modeling/table-figure/tab8-turning-point-comparison.html")

