################################################################################
# Table 2. Variable definitions and model covariates
################################################################################

library(tibble)
library(gt)

# Define the variable table
var_defs <- tribble(
  ~Variable,                ~Definition,                                              ~Type,          ~`Used in`,
  # --- Outcome variables ---
  "count",                  "Total daily trip count",                                  "Count (0+)",   "Models 1-2",
  "count_bin",              "Trip-count bin: 0, 1-2, 3-4, 5-6, 7+",                  "Ordinal (5 levels)", "Model 2",
  "work_day",               "Binary: made any work-related trip on survey day",        "Binary (0/1)", "Model 3a",
  "non_work_trip",          "Non-work daily trip count (total minus work trips)",       "Count (0+)",   "Model 3b",
  # --- Treatment variable ---
  "telework_hour",          "Hours spent teleworking on survey day (0-8)",              "Continuous",   "All models",
  "telework_hour^2",        "Quadratic term for telework hours",                       "Continuous",   "All models",
  # --- Temporal controls ---
  "travel_dow",             "Day of week of the survey travel day",                    "Categorical",  "All models",
  "travel_date_season",     "Season: Spring, Summer, Fall, Winter",                    "Categorical",  "All models",
  "year",                   "Survey wave: 2019, 2021, 2023",                           "Categorical",  "All models",
  # --- Demographic controls ---
  "age",                    "Age group of respondent",                                 "Categorical",  "All models",
  "gender",                 "Gender: Male, Female, Prefer not to answer",              "Categorical",  "All models",
  "education",              "Education: Less than college, College, Graduate, Missing", "Categorical",  "All models",
  # --- Household controls ---
  "num_people",             "Number of people in household",                           "Continuous",   "All models",
  "num_vehicles",           "Number of vehicles in household",                         "Continuous",   "All models",
  "res_type",               "Residence type (single-family, multi-unit, etc.)",        "Categorical",  "All models",
  "income_broad",           "Household income bracket",                                "Categorical",  "All models",
  # --- Spatial controls ---
  "thrive_community_type",  "Community type: Urban, Suburban, Rural, Outside",         "Categorical",  "All models",
  # --- Weights & clustering ---
  "day_weight",             "ACS-rescaled survey weight (divided by 3 for pooling)",   "Weight",       "All models",
  "person_id",              "Person identifier for cluster-robust SEs / survey design","Cluster ID",   "All models"
)

# Create gt table
tab2 <- gt(var_defs) %>%
  tab_header(
    title = "Variable definitions and model covariates"
  ) %>%
  tab_row_group(label = "Outcome Variables",     rows = 1:4) %>%
  tab_row_group(label = "Treatment Variable",    rows = 5:6) %>%
  tab_row_group(label = "Temporal Controls",     rows = 7:9) %>%
  tab_row_group(label = "Demographic Controls",  rows = 10:12) %>%
  tab_row_group(label = "Household Controls",    rows = 13:16) %>%
  tab_row_group(label = "Spatial Controls",      rows = 17) %>%
  tab_row_group(label = "Weights & Clustering",  rows = 18:19) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  cols_width(
    Variable   ~ px(180),
    Definition ~ px(350),
    Type       ~ px(130),
    `Used in`  ~ px(100)
  ) %>%
  tab_options(
    table.font.size = px(11),
    heading.title.font.size = px(14)
  )

# Save
gt::gtsave(tab2, "modeling/table-figure/tab2-variable-definitions.html")


