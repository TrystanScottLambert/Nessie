
b0 <-  100.
r0 <-  180.
ras <-  c(20., 20., 20., -20., -20., -20., 100., 100., 0., 180.)
decs <-  c(-50, -50, -50, 0., 0., 0., 90, 90, 45, -45)
redshifts <-  c(0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.4, 0.4, 0.2, 0.2)
random_redshifts <- rnorm(20000, 0.2, 0.1)
random_redshifts <- random_redshifts[random_redshifts>0]
cosmo <- create_flat_cosmology(0.7, 0.3)
density_function <-  create_density_function(random_redshifts, 20000, 0.001, cosmo)


test_that("getting group ids works on simple case", {
  result <- get_group_ids(b0, r0, ras, decs, redshifts, density_function, cosmo)
  ans <- c(1, 1, 1, 2, 2, 2, 3, 3, -1, -1)
  expect_equal(result, ans)

})
