#' mrgsolve model for a drug
#'
#' The model has the same structure as the Stan model used for estimation and
#' takes the same parameter names, so the posterior estimates from
#' [hmctdmrest()] can be plugged in via `mrgsolve::idata_set()`.
#'
#' @inheritParams get_sample_data
#' @return an `mrgmod` object
#' @export
get_mrgmod <- function(drug) {
  mrgsolve::mread_cache(drug, mrglib(), quiet = TRUE)
}

# Simulate the posterior predictive concentration-time profile:
# one mrgsolve individual per posterior draw, dosed with the patient's events.
# Returns a data.frame with ID (= draw), time, IPRED.
simulate_posterior <- function(hmctdm, n_draws = 200, end = NULL, delta = 0.1) {
  draws <- hmctdm$draws
  if (nrow(draws) > n_draws) draws <- draws[sample.int(nrow(draws), n_draws), , drop = FALSE]
  draws$ID <- seq_len(nrow(draws))
  # patient covariates (LBW, CrCL, ...) are not sampled: take them from idata
  covariates <- setdiff(names(hmctdm$idata), c("ID", names(draws)))
  draws <- cbind(draws, hmctdm$idata[rep(1, nrow(draws)), covariates, drop = FALSE], row.names = NULL)

  events <- hmctdm$data %>%
    filter(.data$evid != 0) %>%
    select(any_of(c("time", "evid", "amt", "cmt", "ss", "ii", "addl", "rate")))
  # cover the whole regimen: last record or last addl dose, whichever is later
  if (is.null(end)) end <- max(hmctdm$data$time, hmctdm$data$time + hmctdm$data$addl * hmctdm$data$ii) * 1.2

  hmctdm$mrgmod %>%
    mrgsolve::idata_set(draws) %>%
    mrgsolve::ev(mrgsolve::as.ev(events)) %>%
    mrgsolve::mrgsim_df(end = end, delta = delta, obsonly = TRUE, recsort = 3)
}
