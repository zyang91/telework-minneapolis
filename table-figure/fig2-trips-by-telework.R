################################################################################
# Figure 2. Mean total daily trips and trip-count-bin distribution
#            by telework intensity category
################################################################################



library(ggplot2)
library(patchwork)

# --- Panel A: Weighted mean trips by telework intensity ---
mean_trips <- trip_count %>%
  group_by(telework_intensity) %>%
  summarise(
    mean_trips = weighted.mean(count, day_weight, na.rm = TRUE),
    se_trips   = sqrt(
      sum(day_weight * (count - weighted.mean(count, day_weight, na.rm = TRUE))^2,
          na.rm = TRUE) /
        (sum(day_weight, na.rm = TRUE)^2 / sum(day_weight^2, na.rm = TRUE))
    ) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

p_a <- ggplot(mean_trips, aes(x = telework_intensity, y = mean_trips)) +
  geom_col(fill = "#4E79A7", width = 0.6) +
  geom_errorbar(aes(ymin = mean_trips - 1.96 * se_trips,
                     ymax = mean_trips + 1.96 * se_trips),
                width = 0.2) +
  geom_text(aes(label = sprintf("%.2f", mean_trips)),
            vjust = -10, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Weighted mean daily trips",
       subtitle = "(a) Mean total daily trips") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 9))

# --- Panel B: Trip-count bin distribution by telework intensity ---
bin_dist <- trip_count %>%
  group_by(telework_intensity, count_bin) %>%
  summarise(wt = sum(day_weight, na.rm = TRUE), .groups = "drop") %>%
  group_by(telework_intensity) %>%
  mutate(pct = 100 * wt / sum(wt)) %>%
  ungroup()

p_b <- ggplot(bin_dist, aes(x = count_bin, y = pct, fill = telework_intensity)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(
    values = c("None (0 hrs)" = "#4E79A7",
               "Moderate (1-6 hrs)" = "#F28E2B",
               "Full-day (7-8 hrs)" = "#E15759"),
    name = "Telework intensity"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "Trip-count bin", y = "Weighted share (%)",
       subtitle = "(b) Trip-count bin distribution") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

fig2 <- p_a / p_b +
  plot_annotation(
    title = "Mean total daily trips and trip-count-bin distribution by telework intensity",
  )

fig2


ggsave("modeling/table-figure/fig2-trips-by-telework.png",
       fig2, width = 7, height = 9.5, dpi = 300)

