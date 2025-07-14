test_that("it runs", {
  ra_array <- c(120., 120., 50.)
  dec_array <- c(-34., -34., 23.)
  redshift_array <- c(0.2, 0.2, 0.6)
  velocity_errors <- c(50., 50., 50.)
  absolute_magnitudes <- c(-18., -18., -18)

  cosmo <- FlatCosmology$new(h = 0.7, omega_m = 0.3)
  completeness <- rep(0.98, length(ra_array))
  density_function <- function(z) {rep(0.2, length(z))}
  cat = RedshiftCatalog$new(ra_array, dec_array, redshift_array, density_function, cosmo, completeness)

  group <- cat$get_raw_groups(0.06, 18)
  cat$run_fof(0.06, 18)
  group_catalog <- cat$calculate_group_table(absolute_magnitudes, velocity_errors)

  expect_equal(group$galaxy_id, c(1, 2))
  expect_equal(group$group_id, c(1, 1))
  expect_equal(cat$group_ids, c(1, 1, -1))
  expect_equal(group_catalog$ra, 120.) # the fact that this is a single number also tests length=1
  expect_equal(group_catalog$dec, -34.)
  expect_equal(group_catalog$multiplicity, 2)
})

test_that("completeness is automatically set", {
  ra_array <- c(120., 120., 50.)
  dec_array <- c(-34., -34., 23.)
  redshift_array <- c(0.2, 0.2, 0.6)
  cosmo <- FlatCosmology$new(h = 0.7, omega_m = 0.3)
  density_function <- function(z) {rep(0.2, length(z))}
  cat <- RedshiftCatalog$new(ra_array, dec_array, redshift_array, density_function, cosmo)
  expect_equal(cat$completeness, rep(1, length(ra_array)))
})


test_that("getting group ids works on simple case", {
  b0 <-  100.
  r0 <-  180.
  ras <-  c(20., 20., 20., 40, 40, 40, 100., 100., 0., 180.)
  decs <-  c(-50, -50, -50, 0., 0., 0., 90, 90, 45, -45)
  redshifts <-  c(0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.4, 0.4, 0.2, 0.2)
  cosmo <- FlatCosmology$new(h = 0.7, omega_m = 0.3)

  random_redshifts <- rnorm(20000, 0.2, 0.1)
  random_redshifts <- random_redshifts[random_redshifts > 0]
  rho_mean <-  create_density_function(random_redshifts, 20000, 0.001, cosmology = cosmo)

  cat <- RedshiftCatalog$new(ra_array = ras, dec_array = decs, redshift_array = redshifts, density_function = rho_mean, cosmology = cosmo)
  cat$run_fof(b0, r0)

  ans <- c(1, 1, 1, 2, 2, 2, 3, 3, -1, -1)
  expect_equal(cat$group_ids, ans)

})

test_that("comparison to mocks is working", {
  group_ids <-  c(1, 1, 1, 2, 2, 2, 2, 3, 3, -1, -1)
  mock_group_ids <-  c(100, 100, 100, 200, 200, -1, -1, -1, -1, -1, -1)
  ras <- rep(120, 11)
  decs <- rep(45, 11)
  redshifts <- rep(0.2, 11)
  cosmo <- FlatCosmology$new(h = 0.7, omega_m = 0.3)
  rho_mean <- function(val){1}
  cat <- RedshiftCatalog$new(ras, decs, redshifts, rho_mean, cosmo)
  cat$group_ids <- group_ids
  cat$mock_group_ids <- mock_group_ids
  metrics <- cat$compare_to_mock()
  expect_equal(metrics, 0.1111, tolerance=1e-3)
})

