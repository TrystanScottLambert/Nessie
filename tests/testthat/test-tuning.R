test_that("tuning works", {
  skip("This test is currently optional because it takes a while.")
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
  lc_numbers <- c(1, 2)
  redshift_catalogues_galform <- list()
  for (lc in lc_numbers) {
    galform_data <- as.data.frame(arrow::read_parquet("~/Desktop/GAMA_paper_plotter/mocks/galform_gals_for_R.parquet"))
    local_catalogue <- galform_data[galform_data$Volume == lc, ]
    local_catalogue <- local_catalogue[local_catalogue$RA < 150, ]

    red_cat <- RedshiftCatalog$new(local_catalogue$RA, local_catalogue$DEC, local_catalogue$Zspec, g09_rho_mean, cosmo)
    red_cat$mock_group_ids <- local_catalogue$GroupID
    red_cat$completeness <- rep(0.95, length(red_cat$ra_array))

    redshift_catalogues_galform[[length(redshift_catalogues_galform) + 1]] <- red_cat
  }


  result <- tune_group_finder(redshift_catalogues_galform, 5, 0.05, 18, c(0.04, 0.09), c(20, 50))
  expect_equal(as.numeric(result$parm['b0']), 0.07010954, tolerance = 1e-5)
  expect_equal(as.numeric(result$parm['r0']), 38.36418, tolerance = 1e-5)
})
