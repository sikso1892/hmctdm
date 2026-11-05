test_that("Cockcroft-Gault divides by serum creatinine and applies the female factor", {
  # (140 - 50) * 60 / (72 * 1.0) = 75
  expect_equal(Cockroft_Gault_CrCl(AGE = 50, SEX = 0, LBW = 60, SCR = 1.0), 75)
  expect_equal(Cockroft_Gault_CrCl(AGE = 50, SEX = 1, LBW = 60, SCR = 1.0), 75 * 0.85)
  expect_equal(Cockroft_Gault_CrCl(AGE = 50, SEX = 0, LBW = 60, SCR = 0.1), 127.1)
  expect_equal(Cockroft_Gault_CrCl(AGE = 50, SEX = 0, LBW = 60, SCR = 20), 10.3)
  # vectorised
  expect_length(Cockroft_Gault_CrCl(AGE = c(50, 80), SEX = c(0, 1), LBW = c(60, 50), SCR = c(1, 2)), 2)
})

test_that("Peck LBW is capped at body weight and vectorised", {
  expect_equal(Peck_LBW(age = 30, sex = 0, wt = 40, ht = 180), 40)
  lbw <- Peck_LBW(age = c(30, 30), sex = c(0, 1), wt = c(90, 90), ht = c(180, 180))
  expect_length(lbw, 2)
  expect_true(all(lbw < 90) && lbw[1] > lbw[2])
})

test_that("b_Cp and the b_Cp dose recommendation are inverses", {
  Vmax <- 500; Km <- 5; CL_R <- 0.1; F <- 0.9; tau <- 12; target <- 15
  dose <- get_recommend_dose_use_b_cp(target, Vmax = Vmax, CL_R = CL_R, Km = Km, F = F, tau = tau)
  expect_equal(get_b_Cp(dose, F = F, tau = tau, CLp = CL_R, Vmax = Vmax, Km = Km), target)
  # no linear clearance: C = r * Km / (Vmax/24 - r)
  dose0 <- get_recommend_dose_use_b_cp(target, Vmax = Vmax, CL_R = 0, Km = Km, F = F, tau = tau)
  expect_equal(get_b_Cp(dose0, F = F, tau = tau, CLp = 0, Vmax = Vmax, Km = Km), target)
  # salt factor round-trips too
  doseS <- get_recommend_dose_use_b_cp(target, Vmax = Vmax, CL_R = CL_R, Km = Km, F = F, tau = tau, salt_factor = 0.92)
  expect_equal(get_b_Cp(doseS, F = F, tau = tau, CLp = CL_R, Vmax = Vmax, Km = Km, salt_factor = 0.92), target)
  # no steady state when input exceeds Vmax
  expect_error(get_b_Cp(10000, F = 1, tau = 12, CLp = 0, Vmax = Vmax, Km = Km), "no steady state")
  # vectorised over dose (and CLp == 0 mixed in)
  v <- get_b_Cp(c(dose, dose0), F = F, tau = tau, CLp = c(CL_R, 0), Vmax = Vmax, Km = Km)
  expect_equal(v, c(target, target))
  expect_length(get_b_Cp(c(100, 200, 300), F = 1, tau = 12, CLp = 0.1, Vmax = Vmax, Km = Km), 3)
})

test_that("get_recommended_dose scales linearly and returns invisibly", {
  out <- capture.output(d <- get_recommended_dose("AUC", target = 400, current_dose = 500, current_status = 300))
  expect_equal(d, 500 / 300 * 400)
  expect_error(get_recommended_dose("foo", target = 1), "unavailable mode")
  expect_error(get_recommended_dose(NULL, target = 1), "unavailable mode")
  expect_error(get_recommended_dose("Cp", target = 1), "current_dose")
})

test_that("b_Cp remains accurate for near-zero linear clearance", {
  expect_equal(get_b_Cp(100, 1, 24, 1e-16, 500, 5), 1.25, tolerance = 1e-12)
})
