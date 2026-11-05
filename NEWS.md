# hmctdmr 0.1.0

## Bug fixes

* Cockcroft-Gault creatinine clearance divided by 72 and then multiplied by
  serum creatinine; it now divides by `72 * SCR`. Estimated CrCL (and hence
  the renal clearance term for amikacin / vancomycin) changes for every
  patient with `SCR != 1`.
* `LBW` / `CrCL` were computed with `ifelse()` on a scalar condition, so only
  the first row's value was used for every record.
* phenytoin Stan model: the log-normal priors were centred on `log(omega)`
  instead of `log(Prior_*)`.
* phenytoin prior names can now be overridden as `Vmax`, `Km`, `VD_NR`
  (previously the names carried a `Prior_` prefix and never matched).
* mrgsolve models for vancomycin and theophiline were copies of the amikacin
  one-compartment model; they are now the two-compartment IV and
  one-compartment oral models matching the Stan models.
* mrgsolve models used the parameter name `CLCR`, but the estimation output
  provides `CrCL`, so simulations silently used the default CrCL of 30 mL/min.
* Stan initial values for vancomycin / theophiline named amikacin parameters;
  inits are now derived from each drug's prior.
* `Vmax` was dropped from the phenytoin posterior summary (the filter pattern
  `x` matched it).
* `get_recommended_dose(mode = "b_Cp")` used `CL_R * (Km - target)`; the
  steady-state mass balance requires `Km + target`. `get_b_Cp()` now solves
  the matching quadratic (`Vd` cancels and was removed from its arguments).
* `data = NULL` fell back to an undefined variable instead of the sample data.
* `print()` swapped the "Stan model option" / "Stan sample option" labels.

## New features

* `plot()` method: posterior predictive concentration-time profile with a
  90% credible band and the observed concentrations.
* `hmctdmrest()` result carries the posterior draws of the PK parameters in
  `$draws`.
* `get_sample_data()` and `get_default_prior()` are exported.
* Input validation: required columns, at least one observation, unknown
  prior names, unknown drug.
* testthat suite that runs without CmdStan (mrgsolve models, helpers,
  post-processing against a mocked fit).

## Housekeeping

* `NAMESPACE` declares imports (the package no longer depends on dplyr being
  attached); S3 methods are registered.
* Manuscript supplementary files are excluded from the built package
  (`.Rbuildignore`).

## Notes

* Default renal clearance slopes (`CL_SLOPE`) and the published validation
  were obtained with the old Cockcroft-Gault expression. Re-check them
  against the corrected CrCL before relying on the defaults clinically.
* For phenytoin, `get_b_Cp()` / `get_recommended_dose(mode = "b_Cp")` expect
  the patient's weight-scaled `Vmax` (`Norm_Vmax` in the mrgsolve output),
  not the 70 kg reference value reported in `$param`.

## Independent review hardening

* Posterior plots now process explicit doses before grid observations at the
  same time (`recsort = 3`), fixing the spurious zero at the start of an `ss = 1`
  regimen and making explicit and `addl` doses consistent. Grid values at a
  dose time are post-dose; observations before a same-time dose remain troughs.
* User-supplied sampler `init` is preserved instead of overwritten by defaults.
* Nested/vector CmdStan model options are retained as single-row list columns,
  preventing postprocessing errors with `cpp_options` or `include_paths`.
* Posterior `cHat` summaries are matched by numeric record index, rather than
  relying on summary row order.
* Multiple patient IDs are rejected: these Stan models fit one patient's
  entire event history, and previously silently pooled multiple patients.
* `get_b_Cp()` uses a stable quadratic root for near-zero linear clearance.
  Its documentation now states that constant-input steady state is only an
  approximation to the average concentration under intermittent nonlinear dosing.
* All four Stan models now use prefix `array` declarations; postfix declarations
  were removed in Stan 2.33. Source regression checks pass, but these changes
  are **untested by compilation or sampling** because CmdStan/Torsten is not
  installed in the test environment. A Torsten-enabled Stan >= 2.26 is required.
* Local `.omx` runtime state is excluded from source package builds.
