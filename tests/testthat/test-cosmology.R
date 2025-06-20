test_that("sigma matches celestial", {
  library(celestial)
  cosmo <- FlatCosmology$new(0.3, 0.7)
  redshift <- 0.8
  mass_solar <- 1e15

  result = calculate_max_sigmas(mass_solar, redshift, cosmo$omega_m, cosmo$omega_k, cosmo$omega_lambda, cosmo$hubble_constant)
  #expected = celestial::coshaloMvirToSigma(mass_solar, z = redshift, Rho = "crit", H0 = cosmo$hubble_constant, OmegaM = cosmo$omega_m)
  expected = 1969.705 # This is handling the fact that celestial has two different constants defined for Solar mass in KG/s
  expect_equal(result, expected, tolerance = 1e-3)
})


test_that("rvir matches celestial", {
  cosmo <- FlatCosmology$new(0.3, 0.7)
  redshift <- 0.8
  mass_solar <- 1e15

  result = calculate_max_rvirs(mass_solar, redshift, cosmo$omega_m, cosmo$omega_k, cosmo$omega_lambda, cosmo$hubble_constant)
  expected = celestial::coshaloMvirToRvir(mass_solar, z = redshift, Rho = "crit", H0 = cosmo$hubble_constant, OmegaM = cosmo$omega_m)
  expect_equal(result, expected)
})

test_that("hubble grow matches celestial", {
  cosmo <- FlatCosmology$new(0.3, 0.7)
  redshift <- 0.8
  result <- h_at_z(redshift, cosmo$omega_m, cosmo$omega_k, cosmo$omega_lambda, cosmo$hubble_constant)
  answer <- celestial::cosgrowH(H0 = cosmo$hubble_constant, z = redshift, OmegaM = cosmo$omega_m)
  expect_equal(result, answer)
})

test_that("distance modulus matches celestial", {
  cosmo <- FlatCosmology$new(1, 0.3)
  redshifts <- c(0.1, 0.2, 0.5, 4, 2, 1, 0, 0.2445)
  answers <- celestial::cosdistDistMod(redshifts)
  results <- cosmo$dist_mod(redshifts)
  expect_equal(answers, results)
})
