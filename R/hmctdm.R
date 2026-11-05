#' @export
print.hmctdm <- function(x, ...) {
  cat("\nStatus: ", x$status, "\n")
  cat("Drug: ", stringr::str_to_title(x$drug), "\n\n")

  cat("\nStan model option: \n"); print(x$stan_model_option)
  cat("\nStan sample option: \n"); print(x$stan_sample_option)
  cat("\nPrior Information: \n"); print(x$prior)
  cat("\nEstimated PK Parameters: \n"); print(x$param)
  cat("\nEstimated Concentration: \n"); print(x$est)
  cat("\nIndividual Data: \n"); print(x$idata)
  invisible(x)
}

#' Posterior predictive concentration-time plot
#'
#' Simulates the patient's dosing history for a sample of posterior draws and
#' plots the median prediction with a 90% credible band over the observed
#' concentrations.
#'
#' @param x an `hmctdm` object from [hmctdmrest()]
#' @param n_draws number of posterior draws to simulate
#' @param end end of the simulation (h); default 1.2 x last record time
#' @param delta simulation time step (h)
#' @param ... unused
#' @return a `ggplot` object
#' @export
plot.hmctdm <- function(x, n_draws = 200, end = NULL, delta = 0.1, ...) {
  sim <- simulate_posterior(x, n_draws = n_draws, end = end, delta = delta)
  band <- sim %>%
    group_by(.data$time) %>%
    summarise(lo = quantile(.data$IPRED, 0.05), med = median(.data$IPRED), hi = quantile(.data$IPRED, 0.95), .groups = "drop")
  obs <- x$data %>% filter(.data$evid == 0)

  ggplot(band, aes(x = .data$time)) +
    geom_ribbon(aes(ymin = .data$lo, ymax = .data$hi), fill = "steelblue", alpha = 0.25) +
    geom_line(aes(y = .data$med), colour = "steelblue") +
    geom_point(data = obs, aes(x = .data$time, y = .data$DV), size = 2) +
    labs(x = "Time (h)", y = "Concentration (mg/L)",
         title = paste(stringr::str_to_title(x$drug), "posterior prediction"),
         subtitle = "line: posterior median, band: 90% credible interval, points: observed") +
    theme_bw()
}
