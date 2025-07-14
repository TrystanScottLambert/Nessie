test_that("RA validation works", {
  expect_silent(validate(c(0, 180, 359.999), type = "ra"))
  expect_error(validate(c(-1, 12, 12), type = "ra"),
               "RA values must be between 0 and 360")
  expect_error(validate(c(12, 12, 361), type = "ra"),
               "RA values must be between 0 and 360")
})

test_that("Dec validation works", {
  expect_silent(validate(c(-90, 0, 89.9), type = "dec"))
  expect_error(validate(c(-90, 91), type = "dec"),
               "Dec Values must be between -90 and 90")
  expect_error(validate(c(-91, 89), type = "dec"),
               "Dec Values must be between -90 and 90")
})

test_that("Redshift validation works", {
  expect_silent(validate(c(0, 1, 1099), type = "redshift"))
  expect_error(validate(c(-0.1, 0.2), type = "redshift"),
               "Redshifts have negative values!")
  expect_warning(validate(c(0.1, 1500), type = "redshift"),
                 "REDSHIFTS ARE HUGE! Are you sure this is correct?")
})

test_that("Absolute magnitude validation works", {
  expect_silent(validate(c(-20, -30), type = "absolute_mag"))
  expect_warning(validate(c(-3), type = "absolute_mag"),
                 "ABSOLUTE MAGNITUDES LOOK WEIRD. MAKE SURE THEY ARE CORRECT!")
  expect_warning(validate(c(-55), type = "absolute_mag"),
                 "ABSOLUTE MAGNITUDES LOOK WEIRD. MAKE SURE THEY ARE CORRECT!")
})

test_that("Velocity error validation works", {
  expect_silent(validate(c(20, 30), type = "vel_err"))
  expect_warning(validate(c(2001), type = "vel_err"),
                 "Velocity error seems very large. Are units correct?")
  expect_error(validate(c(-55), type = "vel_err"),
                 "Velocity error cannot be negative.")
})

test_that("Completeness validation works", {
  expect_silent(validate(c(0.5, 1), type = "completeness"))
  expect_error(validate(c(1.1, 0,5), type = "completeness"),
               "Completeness must be between 0 and 1")
  expect_error(validate(c(-0.1, 0.5), type = "completeness"),
               "Completeness must be between 0 and 1")
})

test_that("b0 validation works", {
  expect_silent(validate(1.0, type = "b0"))
  expect_error(validate(-0.5, type = "b0"),
               "b0 cannot be negative!")
})

test_that("r0 validation works", {
  expect_silent(validate(5.0, type = "r0"))
  expect_error(validate(-0.1, type = "r0"),
               "R0 cannot be negative!")
})


test_that("validate_scalar passes for valid numeric scalar", {
  expect_silent(validate_scalar(42))
})

test_that("validate_scalar fails for non-numeric input", {
  expect_error(validate_scalar("not a number"), "must be numeric")
  expect_error(validate_scalar(TRUE), "must be numeric")
})

test_that("validate_scalar fails for non-scalar input", {
  expect_error(validate_scalar(c(1, 2)), "must be a scalar")
  expect_error(validate_scalar(numeric(0)), "must be a scalar")
})

test_that("validate_scalar fails for NA, NaN, and Inf", {
  expect_error(validate_scalar(NA_real_), "cannot be NA")
  expect_error(validate_scalar(NaN), "cannot be NaN")
  expect_error(validate_scalar(Inf), "cannot be infinite")
  expect_error(validate_scalar(-Inf), "cannot be infinite")
})

test_that("validate_array passes for valid numeric vector", {
  expect_silent(validate_array(c(1.1, 2.2, 3.3)))
})

test_that("validate_array fails for non-numeric input", {
  expect_error(validate_array(c("a", "b")), "must be numeric")
  expect_error(validate_array(c(TRUE, FALSE)), "must be numeric")
})

test_that("validate_array fails for non-vector input", {
  expect_error(validate_array(matrix(1:4, 2, 2)), "must be a 1D vector")
  expect_error(validate_array(array(1:8, dim = c(2, 2, 2))), "must be a 1D vector")
})

test_that("validate_array fails for NA, NaN, Inf values", {
  expect_error(validate_array(c(1, NA, 2)), "contains NA")
  expect_error(validate_array(c(1, NaN, 2)), "contains NaN")
  expect_error(validate_array(c(1, Inf, 2)), "contains infinite")
  expect_error(validate_array(c(1, -Inf, 2)), "contains infinite")
})

