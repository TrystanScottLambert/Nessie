test_that("comparing to the old code.", {
  infile_master_groups <- '~/Desktop/GAMA_paper_plotter/g09_group_catalog_aaron.csv'
  infile_data <- '~/Desktop/GAMA_paper_plotter/gama_galaxy_catalogs/g09_galaxies.dat'
  infile_randoms <- '~/Desktop/GAMA_paper_plotter/gama_g09_randoms.txt'

  df_data <- as.data.frame(read.csv(infile_data, sep = ' '))
  df_data <- df_data[df_data$Z < 0.5, ]
  df_data <- df_data[df_data$Rpetro < 19.8, ]
  velocity_errors <- rep(50., length(df_data$RA))
  cosmo <- create_flat_cosmology(1., 0.3)

  randoms <- as.data.frame(read.csv(infile_randoms))
  random_z <- randoms$z
  g09_frac_area <- 0.001453924
  cosmo_randoms <- create_flat_cosmology(1., 0.25)
  rho_mean <- create_density_function(random_z, length(random_z)/400, g09_frac_area, cosmo_randoms)

  group_ids <- get_group_ids(0.06, 18., df_data$RA, df_data$DEC, df_data$Z, rho_mean, cosmo)
  df_data['group_ids'] <- group_ids
  df_data['ab_mag'] <- df_data$Rpetro - celestial::cosdistDistMod(df_data$Z)

  df_master <- as.data.frame(read.csv(infile_master_groups))
  result <- generate_group_catalog(df_data$RA, df_data$DEC, df_data$Z, df_data$ab_mag, velocity_errors, group_ids, cosmo)

  expect_equal(result$group_id, df_master$Gnum)
  expect_equal(result$multiplicity, df_master$Mult)
  expect_equal(result$ra, df_master$IterCenRA, tolerance = 1e-7) # one or two values incorrect.
  expect_equal(result$dec, df_master$IterCenDEC, tolerance = 1e-7) # one or two values incorrect.
  expect_equal(result$redshift, df_master$MedianZ, tolerance = 1e-4)
  expect_equal(result$co_dist, df_master$MedianDist, tolerance = 1e-3) # just rounding errors
  expect_equal(result$r50, df_master$Rad50, tolerance = 1e-4) # Radii are completely wrong.
  expect_equal(result$r100, df_master$Rad100, tolerance = 1e-4)  # Radii are completely wrong.
  expect_equal(result$rsigma, df_master$Rad1Sig, tolerance = 1e-4)  # Radii are completely wrong.

})
