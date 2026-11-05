#' Bayesian estimation of individual PK parameters
#'
#' Fits the drug's Stan/Torsten model to the patient's dosing and
#' concentration records with HMC, and returns the posterior summaries
#' together with an mrgsolve model parameterised with the posterior medians.
#'
#' @param drug one of `"amikacin"`, `"vancomycin"`, `"theophiline"`, `"phenytoin"`
#' @param data a `data.frame` in the format of [get_sample_data()]; `NULL` uses
#'   the sample data
#' @param prior named list overriding the default prior, see [get_default_prior()]
#' @param stan_model_option list passed to [cmdstanr::cmdstan_model()]
#' @param stan_sample_option list passed to the model's `$sample()` method
#' @return an object of class `hmctdm` with elements `param` (posterior summary),
#'   `est` (data with posterior median `cHat`), `draws` (posterior draws of the
#'   PK parameters), `idata`, `mrgmod`, `prior`, `stan_result`
#' @seealso [plot.hmctdm()], [get_recommended_dose()]
#' @export
hmctdmrest <- function(drug = NULL, data = NULL, prior = NULL,
                       stan_model_option = NULL, stan_sample_option = NULL) {
  object <- list(status = "init") %>%
    preprocess_valid_container(drug, data, prior, stan_model_option, stan_sample_option)

  object$stan_result <- hmctdmr_stan_container(object)
  if (object$status == "init") object$status <- "Estimated"

  hmctdmr_postprocess(object)
}
