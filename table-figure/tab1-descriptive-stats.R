################################################################################
# Table 1. Descriptive statistics of the analytic sample
################################################################################

library(survey)
library(gtsummary)

# Ensure factor types for gtsummary
trip_count <- trip_count %>%
  mutate(
    num_vehicles      = as.numeric(as.character(num_vehicles)),
    travel_dow        = factor(travel_dow),
    age               = factor(age),
    gender            = factor(gender),
    education         = factor(education),
    res_type          = factor(res_type),
    income_broad      = factor(income_broad),
    thrive_community_type = factor(thrive_community_type),
    travel_date_season = factor(travel_date_season),
    year              = factor(year)
  )

# Survey design
svy_design <- svydesign(ids = ~1, weights = ~day_weight, data = trip_count)

# Build Table 1 with survey-weighted descriptive statistics
table1 <- svy_design %>%
  tbl_svysummary(
    by = telework_intensity,
    include = c(count, telework_hour, travel_dow, age, gender, education,
                num_people, num_vehicles, res_type, income_broad,
                thrive_community_type, travel_date_season, year),
    label = list(
      count                 ~ "Daily Trip Count",
      telework_hour         ~ "Telework Hours",
      travel_dow            ~ "Day of Week",
      age                   ~ "Age Group",
      gender                ~ "Gender",
      education             ~ "Education",
      num_people            ~ "Household Size",
      num_vehicles          ~ "Number of Vehicles",
      res_type              ~ "Residence Type",
      income_broad          ~ "Household Income",
      thrive_community_type ~ "Community Type",
      travel_date_season    ~ "Season",
      year                  ~ "Survey Year"
    ),
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous()  ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    type = list(num_vehicles ~ "continuous"),
    missing = "no"
  ) %>%
  add_overall() %>%
  add_n() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Telework Intensity**") %>%
  bold_labels()

# Save as HTML
gt_table1 <- as_gt(table1)
gt::gtsave(gt_table1, "modeling/table-figure/tab1-descriptive-stats.html")


# Print to console
table1

