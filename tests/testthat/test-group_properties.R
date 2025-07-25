test_that("comparing to the old code.", {
  infile_master_groups <- '~/Desktop/GAMA_paper_plotter/g09_group_catalog_aaron.csv'
  infile_data <- '~/Desktop/GAMA_paper_plotter/gama_galaxy_catalogs/g09_galaxies.dat'
  infile_randoms <- '~/Desktop/GAMA_paper_plotter/gama_g09_randoms.txt'

  df_data <- as.data.frame(read.csv(infile_data, sep = ' '))
  df_data <- df_data[df_data$Z < 0.5, ]
  df_data <- df_data[df_data$Rpetro < 19.8, ]
  velocity_errors <- rep(50., length(df_data$RA))
  cosmo <- FlatCosmology$new(h = 1., omega_m = 0.3)

  randoms <- as.data.frame(read.csv(infile_randoms))
  random_z <- randoms$z
  g09_frac_area <- 0.001453924
  cosmo_randoms <- FlatCosmology$new(h = 1., omega_m = 0.25)
  rho_mean <- create_density_function(random_z, length(random_z)/400, g09_frac_area, cosmo_randoms)

  catalog <- RedshiftCatalog$new(df_data$RA, df_data$DEC, df_data$Z, rho_mean, cosmo)
  catalog$set_completeness()
  catalog$run_fof(0.06, 18.)
  group_ids <- catalog$group_ids
  df_data['group_ids'] <- group_ids
  df_data['ab_mag'] <- df_data$Rpetro - celestial::cosdistDistMod(df_data$Z)

  df_master <- as.data.frame(read.csv(infile_master_groups))

  result <- catalog$calculate_group_table(df_data$ab_mag, velocity_errors)

  expect_equal(result$group_id, df_master$Gnum)
  expect_equal(result$multiplicity, df_master$Mult)
  expect_equal(result$iter_ra, df_master$IterCenRA, tolerance = 1e-7) # one or two values incorrect.
  expect_equal(result$iter_dec, df_master$IterCenDEC, tolerance = 1e-7) # one or two values incorrect.
  expect_equal(result$median_redshift, df_master$MedianZ, tolerance = 1e-4)
  expect_equal(result$co_dist, df_master$MedianDist, tolerance = 1e-3) # just rounding errors
  expect_equal(result$r50, df_master$Rad50, tolerance = 1e-4) # Radii are completely wrong.
  expect_equal(result$r100, df_master$Rad100, tolerance = 1e-4)  # Radii are completely wrong.
  expect_equal(result$rsigma, df_master$Rad1Sig, tolerance = 1e-4)  # Radii are completely wrong.

})
