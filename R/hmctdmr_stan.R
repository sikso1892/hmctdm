hmctdmr_stan_container <- function(object) {
  stan_model <- hmctdmr_stan_model(object)
  stanD <- hmctdmr_stan_stanD(object$data, object$prior)
  stanInit <- hmctdmr_stan_stanInit(object$drug, object$prior, object$stan_sample_option$chains)
  fit <- hmctdmr_stan_sample(object, stan_model, stanD, stanInit)

  list(stan_model = stan_model, stan_train = stanD, stan_init = stanInit, fit = fit)
}

hmctdmr_stan_model <- function(object) {
  if (is.null(object$stan_model_option$stan_file)) {
    object$stan_model_option$stan_file <- file.path(stanLib(), paste0(object$drug, ".stan"))
  }
  dir.create(object$stan_model_option$dir, recursive = TRUE, showWarnings = FALSE)
  do.call(cmdstanr::cmdstan_model, object$stan_model_option)
}

hmctdmr_stan_stanD <- function(data, prior) {
  iObs <- which(data$evid == 0)
  c(with(data, list(
    nt = nrow(data), nObs = length(iObs),
    iObs = as.array(iObs), cObs = as.array(DV[iObs]),
    cmt = as.array(cmt), time = as.array(time),
    amt = as.array(amt), rate = as.array(rate),
    evid = as.array(evid), ii = as.array(ii),
    addl = as.array(addl), ss = as.array(ss),
    WT = as.array(WT), LBW = as.array(LBW),
    CrCL = as.array(CrCL)
  )), prior)
}

# Initial values: structural parameters start at their prior centre.
# ka (theophiline) is fixed data in the Stan model, not a parameter.
# phenytoin samples on the log scale (ln_Vmax, ln_Km, ln_VD_NR).
hmctdmr_stan_stanInit <- function(drug, prior, chains) {
  x <- get_prior(drug)
  keep <- !grepl("_omega$", x$name) & x$name != "ka"
  init <- setNames(prior[x$stan_prior_name[keep]], x$name[keep])
  if (drug == "phenytoin") {
    init <- setNames(lapply(init, log), paste0("ln_", names(init)))
  }
  rep(list(init), chains)
}

hmctdmr_stan_sample <- function(object, stan_model, stanD, stanInit = NULL) {
  opt <- object$stan_sample_option
  opt$data <- stanD
  if (is.null(opt$refresh)) opt$refresh <- 0
  if (is.null(opt$init) && !is.null(stanInit)) opt$init <- stanInit
  if (!is.null(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

  suppressMessages(do.call(stan_model$sample, opt))
}
