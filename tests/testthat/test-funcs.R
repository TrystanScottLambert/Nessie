
test_that("density function simple", {
  # This test is difficult to do. We can at least test that this is working and is reasonable.
  # This would also just be useful to have for bench marks.
  # This isn't showing that the function is correct, only that it is unchanged!
  # It would be nice if someone would do a check.

  set.seed(123)  # for reproducibility
  N = 100000
  cosmo <- FlatCosmology$new(0.7, 0.3)
  fractional_area = 0.001
  points <- rnorm(n = N, mean = 0.2, sd = 0.01)
  rho_mean <-  create_density_function(points, N, fractional_area, cosmo)

  redshifts <- c(0., 0.1, 0.3)
  answers <- c(5.386810e-15, 3.898299e-18, 4.453681e-19)

  result = rho_mean(redshifts)
  expect_equal(result, answers)
})

test_that("Group finding graph method simple.", {
  # wild example with very obvious groups just to test core functionality
  ras <-  c(121.1, 121.2, 121.3, 181.1, 181.2, 0.)
  decs <-  c(-23.1, -23.1, -23.1, 68., 68., -30)
  distances <-  c(20., 20., 20., 100., 100., 200.)
  linking_lengths <-  c(2., 2., 2., 2., 2., 2.)
  b0 = rep(0.0349, length(ras))
  r0 = rep(0.4, length(ras))

  result = .find_groups(ras, decs, distances, b0, r0)
  result <-  result[order(result$galaxy_id), ]
  expect_equal(result$galaxy_id, c(1, 2, 3, 4, 5))
  expect_equal(result$group_id, c(1, 1, 1, 2, 2))
})
