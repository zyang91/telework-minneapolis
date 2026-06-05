# Partial and Full telework as Distinct Travel Behaviors: Evidence from Daily Trip-Making in the Twin Cities

Open-source analysis code for a paper on how telework intensity is associated with daily travel outcomes in the Minneapolis–St. Paul region.

*Note: Some preprocessing steps have been simplified or anonymized to protect data confidentiality. The repository provides the main analytical workflow used for the article.*

## Repository purpose

This repository contains scripts used to:

- build pooled analytic datasets across 2019, 2021, and 2023 waves of the MSP Travel Behavior Inventory (TBI),
- estimate main models (ordered logit, negative binomial, mechanism models),
- run propensity score matching (PSM) robustness analyses, and
- generate publication-ready tables and figures.

## Repository structure

```text
telework-minneapolis/
├── psm/
│   ├── 2019-psm.R
│   ├── 2021-psm.R
│   └── 2023-psm.R
└── table-figure/
    ├── 00-common-data.R
    ├── fig1-telework-distribution.R
    ├── fig2-trips-by-telework.R
    ├── fig3-predicted-bin-probs.R
    ├── fig4-predicted-trips-nb.R
    ├── fig5-work-nonwork-mechanism.R
    ├── fig6-psm-balance-bins.R
    ├── fig7-bootstrap-turning-point.R
    ├── tab1-descriptive-stats.R
    ├── tab2-variable-definitions.R
    ├── tab3-ordered-logit.R
    ├── tab4-nb-estimates.R
    ├── tab5-mechanism-models.R
    ├── tab6-psm-robustness.R
    ├── tab7-model-fit-sensitivity.R
    └── tab8-turning-point-comparison.R
```

## Data requirements

The scripts expect cleaned person-day level TBI files for each wave.  
Input paths are intentionally placeholders (e.g., `"YOUR FILE PATH"`), so you can point to your local secure data location.

Key expected variables include:

- outcomes: `count`, `count_bin`, `work_day`, `non_work_trip`
- treatment: `telework_hour`, `telework_intensity`, `telework_binary`
- controls: `travel_dow`, `travel_date_season`, `age`, `gender`, `education`,
  `num_people`, `num_vehicles`, `res_type`, `income_broad`,
  `thrive_community_type`, `year`
- weights/ids: `day_weight`, `person_id`, `day_id`

## Software and R packages

Scripts rely on base R plus packages used throughout the repo, including:

- `tidyverse`, `lubridate`
- `survey`, `MASS`, `ordinal`
- `MatchIt`, `cobalt`
- `gt`, `gtsummary`
- `patchwork`, `tidyr`

Install as needed, for example:

```r
install.packages(c(
  "tidyverse", "lubridate", "survey", "MASS", "ordinal",
  "MatchIt", "cobalt", "gt", "gtsummary", "patchwork", "tidyr"
))
```

## Recommended execution flow

1. **Set file paths** in `table-figure/00-common-data.R` and each `psm/*.R` script where placeholders appear.
2. **Run shared data setup**:
   - source `table-figure/00-common-data.R`
   - this creates `trip_count`, `dat_psm`, and helper functions used downstream.
3. **Run tables/figures scripts** in `table-figure/` as needed.
4. **Run year-specific PSM scripts** in `psm/` for matching diagnostics and robustness checks.

## Script annotations (quick reference)

- `table-figure/00-common-data.R`: pooled data loading/harmonization and helper functions.
- `table-figure/tab1-*.R` to `tab8-*.R`: model tables, variable definitions, robustness and sensitivity outputs.
- `table-figure/fig1-*.R` to `fig7-*.R`: descriptive and model-based visualization scripts.
- `psm/2019-psm.R`, `psm/2021-psm.R`, `psm/2023-psm.R`: year-specific nearest-neighbor matching workflows and balance evaluation.

## Outputs

Most scripts save HTML tables or PNG figures using `gtsave()` / `ggsave()`.  
If needed, update output directories in each script to match your local project layout before running.
