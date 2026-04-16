################################################################################
# Figure 4. Predicted total daily trips across telework intensity levels
#            from the negative binomial model
################################################################################


library(ggplot2)
library(MASS)

# --- 1. Fit quadratic NB model ---
nb_model_q <- glm.nb(
  count ~ travel_dow + age + gender + education +
    telework_hour + I(telework_hour^2) +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  data = trip_count,
  weights = trip_count$day_weight
)

# --- 2. Sample-average predictions across telework grid ---
tw_grid <- seq(0, 8, by = 0.5)
dat <- trip_count
w   <- dat$day_weight

pred_mean <- sapply(tw_grid, function(tw) {
  nd <- dat
  nd$telework_hour <- tw
  mu <- predict(nb_model_q, newdata = nd, type = "response")
  sum(w * mu, na.rm = TRUE) / sum(w, na.rm = TRUE)
})

# --- 3. Parametric bootstrap CI ---
set.seed(123)
B <- 500
b_hat <- coef(nb_model_q)
V_hat <- vcov(nb_model_q)
b_draws <- MASS::mvrnorm(n = B, mu = b_hat, Sigma = V_hat)

X_from_data <- function(nd) model.matrix(formula(nb_model_q), data = nd)

sim_pred <- matrix(NA, nrow = B, ncol = length(tw_grid))
for (j in seq_along(tw_grid)) {
  nd <- dat
  nd$telework_hour <- tw_grid[j]
  Xj <- X_from_data(nd)
  eta_draws <- Xj %*% t(b_draws)
  mu_draws  <- exp(eta_draws)
  sim_pred[, j] <- colSums(mu_draws * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
}

pred_ci <- data.frame(
  telework_hour = tw_grid,
  pred_mean     = pred_mean,
  ci_low  = apply(sim_pred, 2, quantile, probs = 0.025),
  ci_high = apply(sim_pred, 2, quantile, probs = 0.975)
)

# --- 4. Turning point ---
b1 <- coef(nb_model_q)["telework_hour"]
b2 <- coef(nb_model_q)["I(telework_hour^2)"]
tp <- -b1 / (2 * b2)

# --- 5. Plot ---
fig4 <- ggplot(pred_ci, aes(x = telework_hour, y = pred_mean)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2, fill = "#4E79A7") +
  geom_line(linewidth = 1, color = "#4E79A7") +
  geom_vline(xintercept = tp, linetype = "dashed", color = "grey40") +
  annotate("text", x = tp + 0.3, y = max(pred_ci$pred_mean),
           label = sprintf("h* = %.1f hrs", tp),
           hjust = 0, size = 3.5, color = "grey30") +
  scale_x_continuous(breaks = 0:8) +
  labs(
    x = "Telework hours per day",
    y = "Predicted mean total daily trips",
    title = "Predicted total daily trips across telework intensity levels",
    subtitle = "Negative binomial model with quadratic telework; 95% parametric bootstrap CI"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    panel.grid.minor = element_blank()
  )

fig4

ggsave("modeling/table-figure/fig4-predicted-trips-nb.png",
       fig4, width = 7, height = 5, dpi = 300)

cat("Turning point h*:", round(tp, 2), "hours\n")
cat("\nPredicted trips at key levels:\n")
print(pred_ci[pred_ci$telework_hour %in% c(0, 2, 4, 6, 8), ])
cat("\nFigure 4 saved.\n")
