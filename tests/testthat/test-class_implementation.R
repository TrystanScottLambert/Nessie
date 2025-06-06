test_that("it runs", {
  ra_array <- c(120., 120., 50.)
  dec_array <- c(-34., -34., 23.)
  redshift_array <- c(0.2, 0.2, 0.6)
  velocity_errors <- c(50., 50., 50.)
  absolute_magnitudes <- c(-18., -18., -18)

  cosmo <- create_flat_cosmology(0.7, 0.3)
  completeness <- rep(0.98, length(ra_array))
  density_function <- function(z) {rep(0.2, length(z))}
  cat = RedshiftCatalog$new(ra_array, dec_array, redshift_array, density_function, cosmo, completeness)

  group <- cat$fof(0.06, 18)
  cat$set_group_ids(0.06, 18)
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
  cosmo <- create_flat_cosmology(0.7, 0.3)
  density_function <- function(z) {rep(0.2, length(z))}
  cat <- RedshiftCatalog$new(ra_array, dec_array, redshift_array, density_function, cosmo)
  expect_equal(cat$completeness, rep(1, length(ra_array)))
})
