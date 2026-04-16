################################################################################
# Figure 3. Predicted probabilities of daily trip-count bins
#            at 0, 4, and 8 telework hours
################################################################################

library(ggplot2)
library(survey)
library(tidyr)

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

# --- 3. Build reference profile (modal/median) ---
base <- trip_count_cz %>%
  summarise(
    travel_dow            = mode1(travel_dow),
    age                   = mode1(age),
    gender                = mode1(gender),
    education             = mode1(education),
    num_people            = median(num_people, na.rm = TRUE),
    num_vehicles          = median(num_vehicles, na.rm = TRUE),
    res_type              = mode1(res_type),
    income_broad          = mode1(income_broad),
    thrive_community_type = mode1(thrive_community_type),
    travel_date_season    = mode1(travel_date_season),
    work_day              = mode1(work_day),
    year                  = mode1(year)
  )

# --- 4. Convert telework hours to z-scale ---
telework_to_z <- function(h, dat = trip_count_cz) {
  c0 <- h - mean(dat$telework_hour, na.rm = TRUE)
  (c0 - mean(dat$telework_c, na.rm = TRUE)) / sd(dat$telework_c, na.rm = TRUE)
}

# --- 5. Predict at 0, 4, 8 hours ---
scen_048 <- base %>%
  tidyr::crossing(telework_hour = c(0, 4, 8)) %>%
  mutate(
    telework_z  = sapply(telework_hour, telework_to_z),
    telework_z2 = telework_z^2
  )

pp_048 <- predict(ord_mod_z, newdata = scen_048, type = "probs")
tab_048 <- bind_cols(scen_048["telework_hour"], as.data.frame(pp_048))

# --- 6. Also predict on a fine grid for smooth lines ---
grid <- base %>%
  tidyr::crossing(telework_hour = seq(0, 8, by = 0.25)) %>%
  mutate(
    telework_z  = sapply(telework_hour, telework_to_z),
    telework_z2 = telework_z^2
  )

pp_grid <- predict(ord_mod_z, newdata = grid, type = "probs") %>% as.data.frame()

plotdf <- bind_cols(grid["telework_hour"], pp_grid) %>%
  pivot_longer(-telework_hour, names_to = "bin", values_to = "prob") %>%
  mutate(bin = factor(bin, levels = c("0", "1-2", "3-4", "5-6", "7+")))

# Reshape 0/4/8 predictions for point overlay
pts_df <- tab_048 %>%
  pivot_longer(-telework_hour, names_to = "bin", values_to = "prob") %>%
  mutate(bin = factor(bin, levels = c("0", "1-2", "3-4", "5-6", "7+")))

# --- 7. Plot ---
fig3 <- ggplot(plotdf, aes(x = telework_hour, y = prob, color = bin, linetype = bin)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = pts_df, aes(shape = factor(telework_hour)),
             size = 2.5, show.legend = TRUE) +
  scale_color_manual(
    values = c("0" = "#4E79A7", "1-2" = "#F28E2B", "3-4" = "#59A14F",
               "5-6" = "#E15759", "7+" = "#B07AA1"),
    name = "Trip-count bin"
  ) +
  scale_linetype_manual(
    values = c("0" = "solid", "1-2" = "dashed", "3-4" = "dotdash",
               "5-6" = "longdash", "7+" = "twodash"),
    name = "Trip-count bin"
  ) +
  scale_shape_manual(
    values = c("0" = 16, "4" = 17, "8" = 15),
    name = "Telework hours",
    labels = c("0 hrs", "4 hrs", "8 hrs")
  ) +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Telework hours per day",
    y = "Predicted probability",
    title = "Predicted probabilities of daily trip-count bins at 0, 4, and 8 telework hours"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  guides(color = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 1, nrow = 1),
         shape = guide_legend(order = 2, nrow = 1))

fig3

ggsave("modeling/table-figure/fig3-predicted-bin-probs.png",
       fig3, width = 8, height = 5.5, dpi = 300)

cat("\nPredicted bin probabilities at 0/4/8 hours:\n")
print(tab_048)
cat("\nFigure 3 saved.\n")
