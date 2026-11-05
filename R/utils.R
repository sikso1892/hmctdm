#' @importFrom dplyr %>% filter mutate select any_of group_by summarise
#' @importFrom tibble tibble as_tibble
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point labs theme_bw
#' @importFrom rlang .data
#' @importFrom stats setNames median quantile
#' @importFrom utils read.csv
NULL

#' @noRd
mrglib <- function() system.file("mrgsolve", package = "hmctdmr")

#' @noRd
stanLib <- function() system.file("stan", package = "hmctdmr")

#' @noRd
sampleLib <- function() system.file("sample", package = "hmctdmr")

#' Lean body weight (Peck)
#'
#' @param age years
#' @param sex 0 = male, 1 = female
#' @param wt body weight (kg)
#' @param ht height (cm)
#' @return lean body weight (kg), capped at `wt`
#' @noRd
Peck_LBW <- function(age, sex, wt, ht) {
  ht <- ht / 2.54
  LBW <- ifelse(sex == 0,
    ifelse(age <= 18,
      -59.6035 + (5.2878 * ht) - (0.123939 * ht^2) + (0.00128936 * ht^3),
      -130.736 + (4.064 * ht)),
    ifelse(age <= 18,
      -77.55796 + (6.93728 * ht) - (0.171703 * ht^2) + (0.001726 * ht^3),
      -111.621 + (3.636 * ht))
  ) / 2.204623
  pmin(LBW, wt)
}

#' Creatinine clearance (Cockcroft-Gault)
#'
#' CrCL = (140 - AGE) * LBW / (72 * SCR), times 0.85 for female (SEX = 1),
#' clamped to `[crcl_lower, crcl_upper]` mL/min.
#' @noRd
Cockroft_Gault_CrCl <- function(AGE, SEX, LBW, SCR, crcl_lower = 10.3, crcl_upper = 127.1) {
  crcl <- (140 - AGE) * LBW / (72 * SCR)
  crcl <- ifelse(SEX == 0, crcl, crcl * 0.85)
  pmin(pmax(crcl, crcl_lower), crcl_upper)
}

cv_to_sd <- function(cv) sqrt(log(cv^2 + 1))

#' Recommended dose for a target exposure
#'
#' For `mode = "AUC"` or `"Cp"` the dose is scaled linearly:
#' `current_dose / current_status * target`.
#' For `mode = "b_Cp"` (average steady-state concentration under
#' Michaelis-Menten elimination, approximated by constant input) the dose per interval is
#' `(Vmax/24 + CL_R * (Km + target)) / ((Km + target) * F) * target * tau`.
#'
#' @param mode one of `"AUC"`, `"Cp"`, `"b_Cp"`
#' @param target desired value (AUC, concentration, or average concentration)
#' @param current_dose current dose amount
#' @param current_status value of the target quantity under the current regimen
#' @param ... for `mode = "b_Cp"`: `Vmax` (mg/day, for the patient; for phenytoin
#'   use the weight-scaled `Vmax * (WT/70)^0.6`, i.e. the `Norm_Vmax` capture of the
#'   mrgsolve model), `CL_R` (L/h), `Km` (mg/L), `F` (bioavailability),
#'   `tau` (dosing interval, h), `salt_factor` (default 1)
#' @return the recommended dose (invisibly); a summary is printed
#' @export
get_recommended_dose <- function(mode, target, current_dose = NULL, current_status = NULL, ...) {
  available_mode <- c("Cp", "b_Cp", "AUC")
  if (missing(mode) || is.null(mode) || !isTRUE(mode %in% available_mode)) {
    stop(sprintf('unavailable mode("%s") \n-- available mode list: %s',
                 if (missing(mode) || is.null(mode)) "" else toString(mode), toString(available_mode)))
  }

  msg <- list(mode = mode, target = target)

  if (mode == "b_Cp") {
    args <- list(...)
    msg[names(args)] <- args
    recommended_dose <- do.call(get_recommend_dose_use_b_cp, c(list(target = target), args))
  } else {
    if (is.null(current_dose) || is.null(current_status)) {
      stop("current_dose and current_status are required for mode ", mode)
    }
    msg$current_dose <- current_dose
    msg$current_status <- current_status
    recommended_dose <- (current_dose / current_status) * target
  }

  for (key in names(msg)) cat("-", key, ":", msg[[key]], "\n")
  cat("\n- recommended_dose :", recommended_dose, "\n")

  invisible(recommended_dose)
}

# Steady state: dose rate in = elimination rate out
#   r = CL_R * C + (Vmax/24) * C / (Km + C)
# solved for dose = r * tau / F at C = target.
get_recommend_dose_use_b_cp <- function(target, Vmax, CL_R, Km, F, tau, salt_factor = 1) {
  (Vmax / 24 + CL_R * (Km + target)) / ((Km + target) * F * salt_factor) * target * tau
}

#' Constant-input steady-state concentration under Michaelis-Menten elimination
#'
#' Solves `r = CLp * C + (Vmax/24) * C / (Km + C)` for `C`, where
#' `r = dose * salt_factor * F / tau` is the average input rate (mg/h).
#'
#' @param dose dose amount (mg)
#' @param F bioavailability
#' @param tau dosing interval (h)
#' @param CLp linear (e.g. renal) clearance (L/h); 0 if none
#' @param Vmax maximum elimination rate (mg/day) for the patient (for phenytoin
#'   the weight-scaled `Norm_Vmax`, not the 70 kg reference `Vmax`)
#' @param Km Michaelis constant (mg/L)
#' @param salt_factor salt factor (default 1)
#' @details Exact for constant input; an approximation to the average concentration
#' for intermittent dosing because elimination is nonlinear.
#' @return constant-input steady-state concentration (mg/L)
#' @export
get_b_Cp <- function(dose, F, tau, CLp, Vmax, Km, salt_factor = 1) {
  r <- dose * salt_factor * F / tau
  Vmax_h <- Vmax / 24
  no_ss <- CLp == 0 & r >= Vmax_h
  if (any(no_ss)) {
    stop("dose rate (", toString(signif(r[no_ss], 3)), " mg/h) >= Vmax/24: no steady state without linear clearance")
  }
  # CLp*C^2 + (CLp*Km + Vmax_h - r)*C - r*Km = 0   (vectorised over every argument)
  a <- CLp
  b <- CLp * Km + Vmax_h - r
  c <- -r * Km
  disc <- sqrt(b^2 - 4 * a * c)
  # positive root; rationalised form when b >= 0 to avoid cancellation, and it
  # also covers CLp == 0 (a == 0) where it reduces to r*Km/(Vmax_h - r)
  ifelse(b >= 0, -2 * c / (b + disc), (-b + disc) / (2 * a))
}
