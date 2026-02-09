# ------------------------------------------------------------
# simulate_two_period_panel()
# ------------------------------------------------------------
# Generates a two-period long panel dataset with DID structure:
#   Y = alpha + beta1*X1 + beta2*X2 + time_trend[t] + ATT*(D*post) + eps
#
# Returns a data.frame with:
#   id, time, treated, post, X1, X2, Y
# ------------------------------------------------------------

simulate_two_period_panel <- function(N = 100,
                                      N_T = 50,
                                      time_periods = c(0, 1),
                                      true_ATT = 20,
                                      alpha = 10,
                                      beta1 = 3,
                                      beta2 = 5,
                                      p_X2 = 0.4,
                                      time_trend = c(`0` = 0, `1` = 2),
                                      noise_sd = 2,
                                      seed = 123) {

  stopifnot(N_T < N, length(time_periods) == 2)

  set.seed(seed)

  N_C <- N - N_T

  # Covariates
  X1 <- rnorm(N, mean = 0, sd = 1)
  X2 <- rbinom(N, 1, p_X2)

  # Treatment assignment (first N_T treated)
  treated <- c(rep(1, N_T), rep(0, N_C))

  # Build long panel
  df_long <- expand.grid(
    id = 1:N,
    time = time_periods
  )

  df_long$X1 <- X1[df_long$id]
  df_long$X2 <- X2[df_long$id]
  df_long$treated <- treated[df_long$id]
  df_long$post <- ifelse(df_long$time == max(time_periods), 1, 0)

  # Outcome
  eps <- rnorm(nrow(df_long), mean = 0, sd = noise_sd)

  df_long$Y <- alpha +
    beta1 * df_long$X1 +
    beta2 * df_long$X2 +
    time_trend[as.character(df_long$time)] +
    true_ATT * (df_long$treated * df_long$post) +
    eps

  df_long
}
