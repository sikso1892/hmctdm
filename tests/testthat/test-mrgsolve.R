skip_if_not_installed("mrgsolve")

test_that("mrgsolve models compile and expose the Stan parameter names", {
  params <- list(
    amikacin    = c("LBW", "CrCL", "CL_NR", "VD_NR", "CL_SLOPE"),
    vancomycin  = c("LBW", "CrCL", "CL_NR", "VC_NR", "CL_SLOPE", "k12", "k21"),
    theophiline = c("LBW", "CL_NR", "VC_NR", "ka"),
    phenytoin   = c("WT", "LBW", "CrCL", "Vmax", "Km", "VD_NR")
  )
  for (drug in names(params)) {
    mod <- get_mrgmod(drug)
    expect_true(all(params[[drug]] %in% names(mrgsolve::param(mod))), info = drug)
    expect_true("IPRED" %in% mrgsolve::outvars(mod)$capture, info = drug)
  }
})

test_that("amikacin IV bolus follows one-compartment kinetics", {
  mod <- get_mrgmod("amikacin")
  sim <- mod %>%
    mrgsolve::param(LBW = 50, CrCL = 50) %>%
    mrgsolve::ev(amt = 500) %>%
    mrgsolve::mrgsim_df(end = 12, delta = 1)
  CL <- (0.0417 * 50 + 0.815 * 50) * 0.06
  VD <- 0.27 * 50
  expect_equal(sim$IPRED[sim$time == 6], 500 / VD * exp(-CL / VD * 6), tolerance = 1e-3)
})

test_that("vancomycin peripheral compartment receives drug", {
  sim <- get_mrgmod("vancomycin") %>%
    mrgsolve::ev(amt = 1000, rate = 1000) %>%
    mrgsolve::mrgsim_df(end = 12, delta = 1, Req = "CENT PERIPH IPRED")
  expect_true(sim$PERIPH[sim$time == 2] > 0)
  expect_true(all(diff(sim$IPRED[sim$time >= 2]) < 0))
})

test_that("theophiline oral dose is absorbed from GUT into CENT", {
  sim <- get_mrgmod("theophiline") %>%
    mrgsolve::ev(amt = 200, cmt = 1) %>%
    mrgsolve::mrgsim_df(end = 24, delta = 1, obsonly = TRUE)
  expect_equal(sim$IPRED[sim$time == 0], 0)
  expect_true(which.max(sim$IPRED) > 1)
})

# Fake Stan fit: enough of the cmdstanr fit API for hmctdmr_postprocess().
fake_fit <- function(vars, medians, n = 50) {
  summary <- tibble::tibble(variable = vars, mean = medians, median = medians)
  draws <- as.data.frame(lapply(setNames(medians, vars), function(m) m * exp(rnorm(n, 0, 0.1))))
  list(
    summary = function() summary,
    draws = function(variables, format) draws[, variables, drop = FALSE]
  )
}

test_that("postprocess keeps only PK parameters and feeds covariates into mrgsolve", {
  obj <- list(status = "init", drug = "amikacin") %>%
    preprocess_valid_data(NULL) %>%
    preprocess_valid_env(NULL, list(seed = 12)) %>%
    preprocess_valid_prior(NULL)
  obj$stan_result <- list(fit = fake_fit(
    c("lp__", "CL_NR", "VD_NR", "CL_SLOPE", "theta[1]", "theta[2]", "x[1,1]", "x[1,2]",
      "cHat[1]", "cHat[2]", "cHatObs[1]", "sigmaEPS[1]", "CL", "VD"),
    c(-10, 0.042, 0.33, 0.85, 1.85, 14.7, 500, 400, 0, 28, 28, 4.5, 1.85, 14.7)))

  post <- hmctdmr_postprocess(obj)
  expect_s3_class(post, "hmctdm")
  expect_setequal(post$param$variable, c("CL_NR", "VD_NR", "CL_SLOPE", "CL", "VD"))
  expect_equal(post$est$cHat, c(0, 28))
  expect_equal(colnames(post$draws), post$param$variable)
  expect_equal(post$idata$CrCL, max(obj$data$CrCL))

  # mrgsolve CL must use the patient's CrCL/LBW and the posterior medians
  sim <- post$mrgmod %>% mrgsolve::ev(amt = 500) %>% mrgsolve::mrgsim_df(end = 1, delta = 1)
  expect_equal(sim$CL[1], (0.042 * obj$data$LBW[1] + 0.85 * obj$data$CrCL[1]) * 0.06, tolerance = 1e-6)

  expect_output(print(post), "Estimated PK Parameters")
  expect_true("seed" %in% names(post$stan_sample_option))

  # default simulation window covers addl doses
  obj2 <- obj; obj2$data$addl[1] <- 5; obj2$data$ii[1] <- 12
  sim2 <- simulate_posterior(hmctdmr_postprocess(obj2), n_draws = 3, delta = 1)
  expect_gte(max(sim2$time), 60)

  sim <- simulate_posterior(post, n_draws = 10, end = 5, delta = 1)
  expect_equal(length(unique(sim$ID)), 10)
  # every draw must be simulated with the patient's covariates, not model defaults
  sim <- simulate_posterior(post, n_draws = nrow(post$draws), end = 1, delta = 1)
  expect_equal(sim$VD[sim$time == 1], post$draws$VD_NR * obj$data$LBW[1], tolerance = 1e-6)
  p <- plot(post, n_draws = 10)
  expect_s3_class(p, "ggplot")
})

test_that("theophiline ka prior is passed on to the mrgsolve model", {
  obj <- list(status = "init", drug = "theophiline") %>%
    preprocess_valid_data(NULL) %>%
    preprocess_valid_env(NULL, NULL) %>%
    preprocess_valid_prior(list(ka = 1.5))
  obj$stan_result <- list(fit = fake_fit(c("CL_NR", "VC_NR", "CL", "VC", "cHat[1]", "cHat[2]"),
                                         c(40, 0.5, 1.76, 22, 0, 5)))
  post <- hmctdmr_postprocess(obj)
  expect_equal(post$idata$ka, 1.5)
  sim <- simulate_posterior(post, n_draws = 5, end = 2, delta = 1)
  expect_equal(nrow(sim), 15)
})

test_that("postprocess aligns cHat indices and preserves nested model options", {
  obj <- list(status = "init", drug = "amikacin") %>%
    preprocess_valid_data(get_sample_data("amikacin")) %>%
    preprocess_valid_env(list(cpp_options = list(stan_threads = TRUE, stan_opencl = FALSE),
                              include_paths = c("one", "two", "three")), NULL) %>%
    preprocess_valid_prior(NULL)
  obj$stan_result <- list(fit = fake_fit(c("CL_NR", "VD_NR", "CL_SLOPE", "cHat[2]", "cHat[1]"),
                                       c(0.042, 0.3, 0.8, 28, 0)))
  post <- hmctdmr_postprocess(obj)
  expect_equal(post$est$cHat, c(0, 28))
  expect_equal(nrow(post$stan_model_option), 1)
  expect_equal(post$stan_model_option$cpp_options[[1]], obj$stan_model_option$cpp_options)
  expect_equal(post$stan_model_option$include_paths[[1]], c("one", "two", "three"))
})

test_that("posterior simulation applies explicit and additional steady-state doses consistently", {
  mod <- get_mrgmod("amikacin")
  dat <- get_sample_data("amikacin")
  dat$rate <- 0
  dat$ss[1] <- 1
  dat$ii[1] <- 12
  dat$addl[1] <- 2
  obj <- list(data = dat, mrgmod = mod, draws = data.frame(CL_NR = 0.0417, VD_NR = 0.27, CL_SLOPE = 0.815),
              idata = data.frame(ID = 1, LBW = 50, CrCL = 50))
  sim <- simulate_posterior(obj, end = 24, delta = 12)
  ke <- (0.0417 * 50 + 0.815 * 50) * 0.06 / (0.27 * 50)
  peak <- 500 / (0.27 * 50) / (1 - exp(-ke * 12))
  expect_equal(sim$IPRED, rep(peak, 3), tolerance = 1e-4)
  dat$addl[1] <- 0
  doses <- dat[rep(1, 3), ]; doses$time <- c(0, 12, 24); doses$ss <- c(1, 0, 0)
  obj$data <- doses
  explicit <- simulate_posterior(obj, end = 24, delta = 12)
  expect_equal(explicit$IPRED, sim$IPRED, tolerance = 1e-4)
})

test_that("oral steady state and phenytoin infusion follow the model equations", {
  dat <- get_sample_data("theophiline")
  dat$ss[1] <- 1; dat$addl[1] <- 1
  obj <- list(data = dat, mrgmod = get_mrgmod("theophiline"),
              draws = data.frame(CL_NR = 40, VC_NR = 0.5),
              idata = data.frame(ID = 1, LBW = 50, ka = 0.27))
  sim <- simulate_posterior(obj, end = 12, delta = 1)
  t <- sim$time
  expected <- 200 / 25 * 0.27 / (0.27 - 0.08) *
    (exp(-0.08 * t) / (1 - exp(-0.08 * 12)) - exp(-0.27 * t) / (1 - exp(-0.27 * 12)))
  expect_equal(sim$IPRED, expected, tolerance = 1e-4)
  dat <- preprocess_valid_data(list(drug = "phenytoin"), get_sample_data("phenytoin"))$data
  obj <- list(data = dat, mrgmod = get_mrgmod("phenytoin"),
              draws = data.frame(Vmax = 500, Km = 5, VD_NR = 0.8),
              idata = dat[1, c("ID", "WT", "LBW", "CrCL")])
  sim <- simulate_posterior(obj, end = 2, delta = 0.1)
  # Independent RK4 reference, including the exact infusion stop at 92/50 h.
  amt <- 0; h <- 0.0001; vd <- 0.8 * dat$LBW[1]
  f <- function(a, rate) rate - (0.01 * dat$CrCL[1] * 0.06 +
    (500 * (dat$WT[1]/70)^0.6 / 24) / (5 + a/vd)) * a/vd
  for (i in seq_len(20000)) {
    rate <- if (i <= 18400) 50 else 0
    k1 <- f(amt, rate); k2 <- f(amt + h*k1/2, rate)
    k3 <- f(amt + h*k2/2, rate); k4 <- f(amt + h*k3, rate)
    amt <- amt + h/6 * (k1 + 2*k2 + 2*k3 + k4)
  }
  expect_equal(tail(sim$IPRED, 1), amt/vd, tolerance = 1e-4)
})

test_that("postprocess uses the same covariate records as each Stan model", {
  for (drug in available_drugs) {
    dat <- get_sample_data(drug)
    dat$LBW <- c(40, 60); dat$CrCL <- c(20, 100); dat$WT <- c(50, 80)
    obj <- list(status = "init", drug = drug) %>% preprocess_valid_data(dat) %>%
      preprocess_valid_env(NULL, NULL) %>% preprocess_valid_prior(NULL)
    prior <- get_prior(drug)
    keep <- !grepl("_omega$", prior$name) & prior$name != "ka"
    obj$stan_result <- list(fit = fake_fit(c(prior$name[keep], "cHat[1]", "cHat[2]"),
                                         c(prior$value[keep], 0, 5)))
    post <- hmctdmr_postprocess(obj)
    expect_equal(post$idata$LBW, if (drug == "phenytoin") 40 else 60)
    expect_equal(post$idata$CrCL, if (drug == "phenytoin") 20 else 100)
    if (drug == "phenytoin") {
      expect_true("Vmax" %in% post$param$variable)
      expect_equal(post$idata$WT, 50)
    }
  }
})
