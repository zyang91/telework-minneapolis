################################################################################
# Figure 1. Distribution of telework intensity in the pooled weighted sample
################################################################################


library(ggplot2)
library(survey)

# Survey design for weighted statistics
svy <- svydesign(ids = ~1, weights = ~day_weight, data = trip_count)

# Weighted bar chart of telework intensity categories
bar_data <- trip_count %>%
  group_by(telework_intensity) %>%
  summarise(wt_freq = sum(day_weight, na.rm = TRUE), .groups = "drop") %>%
  mutate(pct = round(100 * wt_freq / sum(wt_freq), 1))

fig1 <- ggplot(bar_data, aes(x = telework_intensity, y = pct / 100)) +
  geom_col(
    fill = "#4E79A7",
    color = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = paste0(pct, "%")),
    vjust = -0.5,
    size = 3.5
  ) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    x = "Telework intensity",
    y = "Proportion",
    title = "Distribution of telework intensity in the pooled weighted sample"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  )

fig1

ggsave("modeling/table-figure/fig1-telework-distribution.png",
       fig1, width = 7, height = 4.5, dpi = 300)


