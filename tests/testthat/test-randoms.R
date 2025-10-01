test_that("gen_random_redshifts basic run works", {
  # make a FlatCosmology object like in Python
  cosmo <- FlatCosmology$new(h = 1.0, omega_matter = 0.3)

  # generate some random redshifts (normal around 0.6, sd 0.3) and keep >0
  set.seed(42)  # for reproducibility
  redshifts <- rnorm(500, mean = 0.6, sd = 0.3)
  redshifts <- redshifts[redshifts > 0]

  # generate magnitudes 15–19
  mags <- runif(length(redshifts), min = 15, max = 19)

  z_lim <- max(redshifts) + 0.1
  maglim <- 19.0

  # call your function
  zs <- gen_random_redshifts(
    redshifts = redshifts,
    mags = mags,
    z_lim = z_lim,
    maglim = maglim,
    n_clone = 100,
    iterations = 2,
    cosmo = cosmo
  )

  # same check as Python: roughly n_clone * length(redshifts)
  expect_equal(
    round(length(zs) / (100 * length(redshifts))),
    1
  )
})
