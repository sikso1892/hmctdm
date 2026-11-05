preprocess_valid_container <- function(object, drug, data, prior, stan_model_option, stan_sample_option) {
  object %>%
    preprocess_valid_drug(drug) %>%
    preprocess_valid_data(data) %>%
    preprocess_valid_env(stan_model_option, stan_sample_option) %>%
    preprocess_valid_prior(prior)
}

available_drugs <- c("amikacin", "vancomycin", "theophiline", "phenytoin")

preprocess_valid_drug <- function(object, drug) {
  if (is.null(drug) || !(drug %in% available_drugs)) {
    stop("Unavailable Drug. available: ", toString(available_drugs))
  }
  object$drug <- drug
  object
}

preprocess_valid_data <- function(object, data) {
  if (is.null(data)) {
    message("data is NULL: using sample data for ", object$drug)
    object$status <- "Sample"
    data <- get_sample_data(object$drug)
  }

  required <- c("ID", "time", "evid", "amt", "cmt", "ss", "ii", "addl", "rate",
                "SEX", "AGE", "WT", "HT", "SCR", "DV")
  missing_cols <- setdiff(required, colnames(data))
  if (length(missing_cols) > 0) {
    stop("data is missing required columns: ", toString(missing_cols))
  }
  if (length(unique(data$ID)) != 1L || anyNA(data$ID)) {
    stop("data must contain exactly one patient ID")
  }
  na_cols <- required[vapply(data[required], anyNA, logical(1))]
  if (length(na_cols) > 0) stop("data has NA in required columns: ", toString(na_cols))
  if (!any(data$evid == 0)) stop("data has no observation records (evid == 0)")

  if (!("LBW" %in% colnames(data))) {
    data <- data %>% mutate(LBW = Peck_LBW(.data$AGE, .data$SEX, .data$WT, .data$HT), .after = "WT")
  }
  if (!("CrCL" %in% colnames(data))) {
    data <- data %>% mutate(CrCL = Cockroft_Gault_CrCl(.data$AGE, .data$SEX, .data$LBW, .data$SCR), .after = "SCR")
  }

  object$data <- data
  object
}

preprocess_valid_prior <- function(object, prior) {
  x <- get_prior(object$drug)

  # accept both the short names and the Stan names returned by get_default_prior()
  is_stan_name <- names(prior) %in% x$stan_prior_name
  names(prior)[is_stan_name] <- x$name[match(names(prior)[is_stan_name], x$stan_prior_name)]

  unknown <- setdiff(names(prior), x$name)
  if (length(unknown) > 0) {
    stop("unknown prior name(s): ", toString(unknown), "\n-- available: ", toString(x$name))
  }

  bad <- names(prior)[!vapply(prior, function(v) is.numeric(v) && length(v) == 1, logical(1))]
  if (length(bad) > 0) stop("prior values must be single numbers: ", toString(bad))

  value <- x$value
  override <- x$name %in% names(prior)
  value[override] <- unlist(prior[x$name[override]])

  object$prior <- setNames(as.list(value), x$stan_prior_name)
  object
}

preprocess_valid_env <- function(object, stan_model_option, stan_sample_option) {
  stan_model_option <- as.list(stan_model_option)
  stan_sample_option <- as.list(stan_sample_option)

  set_default <- function(opt, key, value) {
    if (is.null(opt[[key]])) opt[[key]] <- value
    opt
  }

  stan_model_option <- stan_model_option %>%
    set_default("dir", file.path(getwd(), "stanc")) %>%
    set_default("compile", TRUE)
  stan_sample_option <- stan_sample_option %>%
    set_default("chains", 4) %>%
    set_default("thin", 1)
  stan_sample_option <- stan_sample_option %>%
    set_default("iter_warmup", 2500 * stan_sample_option$thin) %>%
    set_default("iter_sampling", 2500 * stan_sample_option$thin) %>%
    set_default("output_dir", file.path(getwd(), "stan_output")) %>%
    set_default("adapt_delta", 0.95)
  object$stan_model_option <- stan_model_option
  object$stan_sample_option <- stan_sample_option
  object
}

#' Sample data set for a drug
#'
#' @param drug one of `"amikacin"`, `"vancomycin"`, `"theophiline"`, `"phenytoin"`
#' @return a `data.frame` in the format expected by [hmctdmrest()]
#' @export
get_sample_data <- function(drug) {
  read.csv(file.path(sampleLib(), paste0(drug, ".csv")))
}

# name:            user-facing prior name (used to override via `prior = list(...)`)
# stan_prior_name: name of the matching `data` entry in the Stan model
# value:           default
get_prior <- function(drug) {
  switch(drug,
    "amikacin" = list(
      name = c("CL_NR", "VD_NR", "CL_SLOPE",
               "CL_NR_omega", "VD_NR_omega", "CL_SLOPE_omega"),
      stan_prior_name = c("Prior_CL_NR", "Prior_VD_NR", "Prior_CL_SLOPE",
                          "Prior_CL_NR_omega", "Prior_VD_NR_omega", "Prior_CL_SLOPE_omega"),
      value = c(0.0417, 0.27, 0.815,
                cv_to_sd(0.25), cv_to_sd(0.3), cv_to_sd(0.4))
    ),
    "vancomycin" = list(
      name = c("CL_NR", "VC_NR", "CL_SLOPE", "k12", "k21",
               "CL_NR_omega", "VC_NR_omega", "CL_SLOPE_omega", "k12_omega", "k21_omega"),
      stan_prior_name = c("Prior_CL_NR", "Prior_VC_NR", "Prior_CL_SLOPE", "Prior_k12", "Prior_k21",
                          "Prior_CL_NR_omega", "Prior_VC_NR_omega", "Prior_CL_SLOPE_omega",
                          "Prior_k12_omega", "Prior_k21_omega"),
      value = c(0.05, 0.21, 0.75, 1.12, 0.48,
                cv_to_sd(0.15), cv_to_sd(0.2), cv_to_sd(0.33), cv_to_sd(0.25), cv_to_sd(0.25))
    ),
    "theophiline" = list(
      name = c("ka", "CL_NR", "VC_NR",
               "CL_NR_omega", "VC_NR_omega"),
      stan_prior_name = c("Prior_ka", "Prior_CL_NR", "Prior_VC_NR",
                          "Prior_CL_NR_omega", "Prior_VC_NR_omega"),
      value = c(0.27, 40, 0.5,
                cv_to_sd(0.15), cv_to_sd(0.2))
    ),
    "phenytoin" = list(
      name = c("Vmax", "Km", "VD_NR",
               "Vmax_omega", "Km_omega", "VD_NR_omega"),
      stan_prior_name = c("Prior_Vmax", "Prior_Km", "Prior_VD_NR",
                          "Prior_Vmax_omega", "Prior_Km_omega", "Prior_VD_NR_omega"),
      value = c(500, 5.0, 0.8,
                cv_to_sd(0.2), cv_to_sd(0.2), cv_to_sd(0.2))
    ),
    stop("Unavailable Drug: ", drug)
  )
}

#' Default prior for a drug
#'
#' @inheritParams get_sample_data
#' @return a named list of prior values as passed to the Stan model
#' @export
get_default_prior <- function(drug) {
  x <- get_prior(drug)
  setNames(as.list(x$value), x$stan_prior_name)
}
