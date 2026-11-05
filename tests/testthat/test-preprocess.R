test_that("prior override replaces only the named entries", {
  obj <- preprocess_valid_prior(list(drug = "amikacin"), list(CL_NR = 0.05))
  expect_equal(obj$prior$Prior_CL_NR, 0.05)
  expect_equal(obj$prior$Prior_VD_NR, 0.27)
  expect_equal(names(obj$prior), names(get_default_prior("amikacin")))
  expect_error(preprocess_valid_prior(list(drug = "amikacin"), list(bogus = 1)), "unknown prior")
  expect_error(preprocess_valid_prior(list(drug = "amikacin"), list(CL_NR = c(0.05, 0.06), VD_NR = 0.3)), "single numbers")
  expect_error(preprocess_valid_prior(list(drug = "amikacin"), list(CL_NR = "a")), "single numbers")
})

test_that("prior can be overridden with the names returned by get_default_prior()", {
  p <- get_default_prior("amikacin"); p$Prior_CL_NR <- 0.05
  obj <- preprocess_valid_prior(list(drug = "amikacin"), p)
  expect_equal(obj$prior$Prior_CL_NR, 0.05)
})

test_that("phenytoin prior can be overridden by short name", {
  obj <- preprocess_valid_prior(list(drug = "phenytoin"), list(Vmax = 400))
  expect_equal(obj$prior$Prior_Vmax, 400)
})

test_that("sample data gets LBW and CrCL added, NULL data falls back to sample", {
  for (drug in c("amikacin", "vancomycin", "theophiline", "phenytoin")) {
    obj <- preprocess_valid_data(list(drug = drug, status = "init"), NULL)
    expect_equal(obj$status, "Sample")
    expect_true(all(c("LBW", "CrCL") %in% colnames(obj$data)))
    expect_true(all(obj$data$LBW <= obj$data$WT))
  }
  expect_error(preprocess_valid_data(list(drug = "amikacin"), data.frame(ID = 1)), "missing required")
  d <- get_sample_data("amikacin"); d$evid[2] <- NA
  expect_error(preprocess_valid_data(list(drug = "amikacin"), d), "NA in required columns: evid")
})

test_that("stan inits match the parameters block of each model", {
  expected <- list(
    amikacin    = c("CL_NR", "VD_NR", "CL_SLOPE"),
    vancomycin  = c("CL_NR", "VC_NR", "CL_SLOPE", "k12", "k21"),
    theophiline = c("CL_NR", "VC_NR"),
    phenytoin   = c("ln_Vmax", "ln_Km", "ln_VD_NR")
  )
  for (drug in names(expected)) {
    init <- hmctdmr_stan_stanInit(drug, get_default_prior(drug), chains = 2)
    expect_length(init, 2)
    expect_setequal(names(init[[1]]), expected[[drug]])
    stan <- readLines(system.file("stan", paste0(drug, ".stan"), package = "hmctdmr"))
    for (p in expected[[drug]]) expect_true(any(grepl(paste0("real(<lower=0>)? ", p, ";"), stan)), info = paste(drug, p))
  }
  expect_equal(hmctdmr_stan_stanInit("phenytoin", get_default_prior("phenytoin"), 1)[[1]]$ln_Vmax, log(500))
})

test_that("default sampler options are filled in without touching the file system", {
  before <- list.files(getwd())
  obj <- preprocess_valid_env(list(), NULL, list(chains = 2, seed = 7))
  expect_equal(list.files(getwd()), before)
  expect_equal(obj$stan_sample_option$seed, 7)
  expect_equal(obj$stan_sample_option$chains, 2)
  expect_equal(obj$stan_sample_option$iter_warmup, 2500)
  expect_true(obj$stan_model_option$compile)
})

test_that("single-patient models reject pooled patient records", {
  dat <- get_sample_data("amikacin"); dat$ID[2] <- 2
  expect_error(preprocess_valid_data(list(drug = "amikacin"), dat), "one patient ID")
})

test_that("sampler respects caller initial values and forwards options", {
  sampler <- list(sample = function(...) list(...))
  obj <- list(stan_sample_option = list(init = 0, seed = 12, parallel_chains = 2))
  args <- hmctdmr_stan_sample(obj, sampler, list(nt = 2), list(list(CL_NR = 1)))
  expect_identical(args$init, 0)
  expect_equal(args$seed, 12)
  expect_equal(args$parallel_chains, 2)
  expect_equal(args$data, list(nt = 2))
  obj$stan_sample_option$init <- NULL
  expect_equal(hmctdmr_stan_sample(obj, sampler, list(), list(list(CL_NR = 1)))$init,
               list(list(CL_NR = 1)))
})

test_that("Stan sources use supported array declarations", {
  for (drug in available_drugs) {
    stan <- readLines(file.path(stanLib(), paste0(drug, ".stan")))
    expect_false(any(grepl("(real|int)\\[\\]", stan)))
    expect_false(any(grepl("(real|int)(<[^>]+>)? [[:alnum:]_]+\\[[^]]+\\];", stan)))
    expect_true(any(grepl("array\\[nt\\] real time;", stan)))
  }
})
