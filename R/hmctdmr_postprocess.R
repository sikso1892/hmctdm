# Stan variables that are not PK parameters (derived/internal quantities).
not_pk_pattern <- "^(lp__|theta|x\\[|cHat|sigmaEPS|parms|ln_)"

hmctdmr_postprocess <- function(hmctdm) {
  fit <- hmctdm$stan_result$fit
  summary <- fit$summary()
  pk_parm <- summary %>% filter(!stringr::str_detect(.data$variable, not_pk_pattern))
  cHat <- summary %>% filter(stringr::str_detect(.data$variable, "^cHat\\["))

  cHat <- cHat[match(paste0("cHat[", seq_len(nrow(hmctdm$data)), "]"), cHat$variable), ]

  post <- list(status = hmctdm$status, drug = hmctdm$drug)

  option_tibble <- function(name, opt) {
    as_tibble(c(list(option = name),
                lapply(opt, function(x) if (length(x) == 1 && is.atomic(x)) x else list(x))))
  }
  post$stan_model_option <- option_tibble("stan_model_option", hmctdm$stan_model_option)
  post$stan_sample_option <- option_tibble("stan_sample_option", hmctdm$stan_sample_option)
  post$prior <- tibble(parameter = names(hmctdm$prior), value = unlist(hmctdm$prior))
  post$param <- pk_parm
  post$data <- hmctdm$data
  post$est <- hmctdm$data
  post$est$cHat <- cHat$median
  post$draws <- as.data.frame(fit$draws(variables = pk_parm$variable, format = "df"))
  post$draws <- post$draws[, pk_parm$variable, drop = FALSE]

  covariates <- hmctdm$data[1, intersect(c("ID", "SEX", "AGE", "WT", "LBW", "HT", "SCR", "CrCL"),
                                         colnames(hmctdm$data)), drop = FALSE]
  # amikacin/vancomycin/theophiline Stan models use max(LBW), max(CrCL); phenytoin uses the first record
  if (hmctdm$drug != "phenytoin") {
    covariates$LBW <- max(hmctdm$data$LBW)
    covariates$CrCL <- max(hmctdm$data$CrCL)
  }
  # ka is fixed data in the theophiline Stan model, so it is not in the posterior summary
  if (hmctdm$drug == "theophiline") covariates$ka <- hmctdm$prior$Prior_ka
  post$idata <- cbind(covariates, as.list(setNames(pk_parm$median, pk_parm$variable)), row.names = NULL)

  post$mrgmod <- get_mrgmod(hmctdm$drug) %>%
    mrgsolve::idata_set(post$idata) %>%
    mrgsolve::carry_out("cmt, ii, addl, rate, amt, evid, ss, SEX, AGE, WT, HT, LBW, SCR, CrCL")

  post$stan_result <- hmctdm$stan_result

  class(post) <- "hmctdm"
  post
}
