
library(arrow)
library(data.table)
library(celestial)
infile_randoms = '~/Desktop/GAMA_paper_plotter/gama_g09_randoms.txt'
infile_data = '~/Desktop/GAMA_paper_plotter/gama_galaxy_catalogs/g09_galaxies.dat'
b0 = 0.06
r0 = 18.
cosmo_randoms  <-  create_flat_cosmology(1.0, 0.25)
cosmo <-  create_flat_cosmology(1., 0.3)
g09_frac_area <- 0.001453924
Mmax <- 1e15

randoms = as.data.table(read.csv(infile_randoms))
random_z = randoms$z

data = as.data.table(read.csv(infile_data, sep = ' '))
data = data[data$Z < 0.5,]
data = data[data$Rpetro < 19.8,]
small_data = data[1:100, ]

dists = comoving_distances_at_z(data$Z, cosmo$Om0, cosmo$OmK, cosmo$OmL, cosmo$H0)

rho_mean = create_density_function(random_z, length(random_z)/400, g09_frac_area, cosmo_randoms)

test_that("can do gama quickly and equally", {

  start.now = Sys.time()
  groups = fof(b0, r0, data$RA, data$DEC, data$Z, rho_mean, cosmology=cosmo)
  end.now = Sys.time()
  print(paste("current time: ", end.now - start.now))

  start.now = Sys.time()
  groups_classic = .fof_aaron(b0, r0, data$RA, data$DEC, data$Z, rho_mean, cosmology=cosmo)
  end.now = Sys.time()
  print(paste("classic time: ", end.now - start.now))

  sorted_groups = groups[order(groups$galaxy_id), ]
  sorted_classic_groups = groups_classic[order(groups_classic$galaxy_id), ]

  expect_equal(sorted_groups$galaxy_id, sorted_classic_groups$galaxy_id)
  expect_equal(sorted_groups$group_id, sorted_classic_groups$group_id)

})

test_that("Benchmarks", {
  # Calculating the plane-of-sky linking lengths
  b0 = 0.06
  r0 = 18.
  Mmax= 1e15
  cosmology = cosmo
  co_dists = dists
  redshifts = data$Z
  density_function = rho_mean

  linking_lengths = density_function(data$Z)^(-1./3)
  gal_rad = b0 * linking_lengths
  max_on_sky_radius = celestial::coshaloMvirToRvir(Mmax, z = data$Z, Rho="crit", cosmology$H0)
  too_wide =  gal_rad > max_on_sky_radius
  gal_rad[too_wide] = max_on_sky_radius[too_wide]
  linking_lengths_pos = gal_rad/(cosmology$h * co_dists)

  # Calculating the line-of-sight linking lengths
  R = r0 * (1 + redshifts)/(sqrt(cosmology$Om0 * (1+redshifts)^3 + cosmology$OmL))
  linking_lengths_los = (gal_rad * R)/cosmology$h
  max_los_distances = celestial::coshaloMvirToSigma(Mmax, z= redshifts, Rho="crit", cosmology$H0) * (1+redshifts) / celestial::cosgrowH(H0=cosmology$H0, z=redshifts)
  too_far = linking_lengths_los > max_los_distances
  linking_lengths_los[too_far] = max_los_distances[too_far]

  start.raw_links = Sys.time()
  fof_links_fast(data$RA, data$DEC, co_dists,  linking_lengths_pos, linking_lengths_los)
  end.raw_links = Sys.time()





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

  start.with_properties = Sys.time()
  group_ids <- get_group_ids(0.06, 18., df_data$RA, df_data$DEC, df_data$Z, rho_mean, cosmo)
  generate_group_catalog(df_data$RA, df_data$DEC, df_data$Z, df_data$ab_mag, velocity_errors, group_ids, cosmo)
  end.with_properties = Sys.time()

  start.without_properties = Sys.time()
  group_ids <- get_group_ids(0.06, 18., df_data$RA, df_data$DEC, df_data$Z, rho_mean, cosmo)
  end.without_properties = Sys.time()

  print(paste('Raw linking time: ', end.raw_links - start.raw_links))
  print(paste("With properties", end.with_properties - start.with_properties))
  print(paste("Without properties", end.without_properties - start.without_properties))


})
