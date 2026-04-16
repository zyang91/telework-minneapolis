################################################################################
# 00-common-data.R
# Shared data loading, pooling, and harmonization for all figure/table scripts
#
# SOURCE this file at the top of each figure/table script:
#   source("modeling/table-figure/00-common-data.R")
#
# Objects created:
#   trip_count   - pooled & harmonized person-day data (main sample, 3 years)
#   dat_psm      - pooled & harmonized PSM matched data (3 years)
################################################################################

library(tidyverse)
library(lubridate)
options(scipen = 999)

################################################################################
# 1. MAIN SAMPLE: LOAD, POOL, HARMONIZE
################################################################################

trip_count2019 <- read.csv("YOUR FILE PATH FOR 2019 MSP TBI DATA")
trip_count2021 <- read.csv("YOUR FILE PATH FOR 2021 MSP TBI DATA")
trip_count2023 <- read.csv("YOUR FILE PATH FOR 2023 MSP TBI DATA")

# Coerce num_people / num_vehicles to character before binding so
# parse_number() can handle all years uniformly (some years store
# these as integer, others as character).
coerce_hh <- function(df) {
  df %>% mutate(
    num_people   = as.character(num_people),
    num_vehicles = as.character(num_vehicles)
  )
}

trip_count <- bind_rows(
  coerce_hh(trip_count2019) %>% mutate(year = 2019),
  coerce_hh(trip_count2021) %>% mutate(year = 2021),
  coerce_hh(trip_count2023) %>% mutate(year = 2023)
)

# ACS-style weight rescaling for pooled waves
trip_count <- trip_count %>%
  mutate(day_weight = day_weight / 3,
         year = as.factor(year))

# Parse numeric household variables
trip_count <- trip_count %>%
  mutate(
    num_people   = parse_number(num_people),
    num_vehicles = parse_number(num_vehicles)
  )

#  Merge category

################################################################################
# 2. PSM SAMPLE: LOAD, POOL, HARMONIZE
################################################################################



################################################################################
# 3. HELPER: cluster-robust vcov for glm.nb
################################################################################

cluster_vcov_fn <- function(model, cluster) {
  cluster <- as.factor(cluster)
  M <- nlevels(cluster); N <- nobs(model); K <- model$rank
  dfc <- (M / (M - 1)) * ((N - 1) / (N - K))
  U   <- sandwich::estfun(model)
  Uc  <- apply(U, 2, function(x) tapply(x, cluster, sum))
  meat <- crossprod(Uc) / N
  dfc * sandwich::sandwich(model, meat = meat)
}

################################################################################
# 4. HELPER: modal value
################################################################################

mode1 <- function(x) names(sort(table(x), decreasing = TRUE))[1]
