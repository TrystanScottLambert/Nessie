

test_that("can do gama quickly", {
  library(arrow)
  library(data.table)
  library(celestial)
  infile_randoms = '~/Desktop/GAMA_paper_plotter/gama_g09_randoms.txt'
  infile_data = '~/Desktop/GAMA_paper_plotter/gama_galaxy_catalogs/g09_galaxies.dat'
  b0 = 0.06
  r0 = 18.
  cosmo_randoms  <-  create_flat_cosmology(1.0, 0.25)
  cosmo <-  create_flat_cosmology(1., 0.3)
  g09_frac_area = 0.001453924
  Mmax = 1e15

  randoms = as.data.table(read.csv(infile_randoms))
  random_z = randoms$z

  data = as.data.table(read.csv(infile_data, sep = ' '))
  data = data[data$Z < 0.5,]
  data = data[data$Rpetro < 19.8,]
  dists = comoving_distances_at_z(data$Z, cosmo$Om0, cosmo$OmK, cosmo$OmL, cosmo$H0)

  rho_mean = create_density_function(random_z, length(random_z)/400, g09_frac_area, cosmo_randoms)
  groups = fof(b0, r0, data$RA, data$DEC, data$Z, rho_mean, cosmo)



})
