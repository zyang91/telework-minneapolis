################################################################################
# Figure 5. Predicted probability of making any work-related trip and
#            predicted non-work trips across telework intensity levels
################################################################################

library(ggplot2)
library(MASS)
library(survey)
library(patchwork)

# --- 1. Prepare work_day as binary ---
trip_count <- trip_count %>%
  mutate(
    work_day_bin = ifelse(work_day == "yes", 1, 0),
    telework_c   = telework_hour - mean(telework_hour, na.rm = TRUE),
    telework_c2  = telework_c^2
  )

# --- 2A. Work-travel logistic regression ---
des <- svydesign(ids = ~person_id, weights = ~day_weight, data = trip_count)

m_work <- svyglm(
  work_day_bin ~ travel_dow + age + gender + education +
    telework_c + telework_c2 +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  design = des,
  family = quasibinomial()
)

# Predictions on a grid
tw_grid <- seq(0, 8, by = 0.5)
tw_mean <- mean(trip_count$telework_hour, na.rm = TRUE)

nd_work <- data.frame(
  telework_hour = tw_grid,
  telework_c    = tw_grid - tw_mean,
  telework_c2   = (tw_grid - tw_mean)^2,
  num_people    = mean(trip_count$num_people, na.rm = TRUE),
  num_vehicles  = mean(trip_count$num_vehicles, na.rm = TRUE)
)
for (v in names(m_work$xlevels)) {
  nd_work[[v]] <- factor(m_work$xlevels[[v]][1], levels = m_work$xlevels[[v]])
}

tt <- delete.response(terms(m_work))
X_work <- model.matrix(tt, data = nd_work, xlev = m_work$xlevels)
b_work <- coef(m_work)
V_work <- vcov(m_work)

eta_work <- as.numeric(X_work %*% b_work)
se_work  <- sqrt(diag(X_work %*% V_work %*% t(X_work)))

pred_work <- data.frame(
  telework_hour = tw_grid,
  p_hat = plogis(eta_work),
  p_lo  = plogis(eta_work - 1.96 * se_work),
  p_hi  = plogis(eta_work + 1.96 * se_work)
)

# --- 2B. Non-work NB model ---
trip_count$non_work_trip[is.na(trip_count$non_work_trip)] <- 0

nb_nonwork <- glm.nb(
  non_work_trip ~ travel_dow + age + gender + education +
    telework_hour + I(telework_hour^2) +
    num_people + num_vehicles + res_type + income_broad +
    thrive_community_type + travel_date_season + year,
  data = trip_count,
  weights = trip_count$day_weight
)

dat <- trip_count
w   <- dat$day_weight

pred_nw_mean <- sapply(tw_grid, function(tw) {
  nd <- dat
  nd$telework_hour <- tw
  mu <- predict(nb_nonwork, newdata = nd, type = "response")
  sum(w * mu, na.rm = TRUE) / sum(w, na.rm = TRUE)
})

# Bootstrap CI for non-work
set.seed(456)
B <- 500
b_hat_nw <- coef(nb_nonwork)
V_hat_nw <- vcov(nb_nonwork)
b_draws_nw <- MASS::mvrnorm(n = B, mu = b_hat_nw, Sigma = V_hat_nw)
X_nw_fn <- function(nd) model.matrix(formula(nb_nonwork), data = nd)

sim_nw <- matrix(NA, nrow = B, ncol = length(tw_grid))
for (j in seq_along(tw_grid)) {
  nd <- dat
  nd$telework_hour <- tw_grid[j]
  Xj <- X_nw_fn(nd)
  eta_draws <- Xj %*% t(b_draws_nw)
  mu_draws  <- exp(eta_draws)
  sim_nw[, j] <- colSums(mu_draws * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
}

pred_nonwork <- data.frame(
  telework_hour = tw_grid,
  pred_mean = pred_nw_mean,
  ci_low  = apply(sim_nw, 2, quantile, probs = 0.025),
  ci_high = apply(sim_nw, 2, quantile, probs = 0.975)
)

# --- 3. Panel A: Work-travel probability ---
p_a <- ggplot(pred_work, aes(x = telework_hour, y = p_hat)) +
  geom_ribbon(aes(ymin = p_lo, ymax = p_hi), alpha = 0.2, fill = "#E15759") +
  geom_line(linewidth = 1, color = "#E15759") +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Telework hours per day",
       y = "Pr(any work-related trip)",
       subtitle = "(a) Predicted probability of making any work-related trip") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

# --- 4. Panel B: Non-work predicted trips ---
p_b <- ggplot(pred_nonwork, aes(x = telework_hour, y = pred_mean)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2, fill = "#59A14F") +
  geom_line(linewidth = 1, color = "#59A14F") +
  scale_x_continuous(breaks = 0:8) +
  labs(x = "Telework hours per day",
       y = "Predicted mean non-work trips",
       subtitle = "(b) Predicted non-work trips") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

fig5 <- p_a / p_b +
  plot_annotation(
    title = "Predicted probability of work-related trip and non-work trips across telework intensity",
    subtitle = "Survey-weighted logistic (a) and negative binomial (b) models; 95% CI shaded",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey40")
    )
  )

fig5


ggsave("modeling/table-figure/fig5-work-nonwork-mechanism.png",
       fig5, width = 7, height = 8, dpi = 300)

cat("Figure 5 saved.\n")
