test_that("tuning works", {
  library(celestial)
  library(data.table)
  library(arrow)

  # Setting up cosmology
  cosmo <- Nessie::FlatCosmology$new(h = 0.67, omega_matter = 0.3)

  # Calculating rho mean from random catalogues.
  g09_area <- skyarea(c(129, 141), c(-2, 3))
  g12_area <- skyarea(c(174, 186), c(-3, 2))
  g15_area <- skyarea(c(211.5, 223.5), c(-2, 3))
  g23_area <- skyarea(c(339, 351), c(-35, -30))

  g09_randoms <- as.data.table(read.csv("~/Desktop/GAMA_paper_plotter/gama_g09_randoms.txt"))
  g12_randoms <- as.data.table(read.csv("~/Desktop/GAMA_paper_plotter/gama_g12_randoms.txt"))
  g15_randoms <- as.data.table(read.csv("~/Desktop/GAMA_paper_plotter/gama_g15_randoms.txt"))
  g23_randoms <- as.data.table(read.csv("~/Desktop/GAMA_paper_plotter/gama_g23_randoms.txt"))

  g09_rho_mean <- Nessie::create_density_function(g09_randoms$z, length(g09_randoms$z)/400, g09_area["areafrac"], cosmo)
  g12_rho_mean <- Nessie::create_density_function(g12_randoms$z, length(g12_randoms$z)/400, g12_area["areafrac"], cosmo)
  g15_rho_mean <- Nessie::create_density_function(g15_randoms$z, length(g15_randoms$z)/400, g15_area["areafrac"], cosmo)
  g23_rho_mean <- Nessie::create_density_function(g23_randoms$z, length(g23_randoms$z)/400, g23_area["areafrac"], cosmo)
  rho_means <- list(g09 = g09_rho_mean, g12 = g12_rho_mean, g15 = g15_rho_mean, g23 = g23_rho_mean)

  # Setting up Redshift catalogues
  rlim <- 19.65

  calibration_data <- as.data.frame(arrow::read_parquet("~/Desktop/GAMA_paper_plotter/mocks/gama_gals_for_R.parquet"))
  calibration_data <- calibration_data[calibration_data['zobs'] < 0.5, ]


  lightcone_numbers = c(0)#, 3, 4, 5, 7, 8, 9, 10)
  gama_fields = c('g09')#, 'g12', 'g15', 'g23')
  combinations <- expand.grid(lightcone_numbers, gama_fields)
  all_fields <- paste0(combinations$Var1, combinations$Var2)

  redshift_catalogues <- list()
  for (g_field in all_fields) {
    gama_field <- substr(g_field, nchar(g_field) - 2, nchar(g_field))
    local_catalogue <- calibration_data[calibration_data['lightcone_gamafield'] == g_field, ]
    red_cat <- RedshiftCatalog$new(local_catalogue$ra, local_catalogue$dec, local_catalogue$zobs, rho_means[[gama_field]], cosmo)
    red_cat$mock_group_ids <- local_catalogue$GroupID

    # setting this values to -1
    counts <- table(red_cat$mock_group_ids)
    singleton_ids <- names(counts[counts == 1])
    red_cat$mock_group_ids <- ifelse(red_cat$mock_group_ids %in% singleton_ids, -1, red_cat$mock_group_ids)

    redshift_catalogues[[length(redshift_catalogues) + 1]] <- red_cat
  }

  result <- tune_group_finder(redshift_catalogues, 5, 0.05, 18, c(0.03, 0.09), c(1, 100))
  expect_equal(as.numeric(result$parm['b0']), 0.07010954, tolerance = 1e-5)
  expect_equal(as.numeric(result$parm['r0']), 38.36418, tolerance = 1e-5)
})
