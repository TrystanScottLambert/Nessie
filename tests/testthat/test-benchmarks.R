

test_that("can do gama quickly", {
  library(arrow)
  library(data.table)
  library(celestial)
  infile_randoms = '~/Desktop/GAMA_paper_plotter/gama_g09_randoms.txt'
  infile_data = '~/Desktop/GAMA_paper_plotter/gama_galaxy_catalogs/g09_galaxies.dat'
  cosmo = create_flat_cosmology(0.7, 0.3)
  g09_frac_area = 0.001453924
  Mmax = 1e15

  randoms = as.data.table(read.csv(infile_randoms))
  random_z = randoms$z

  data = as.data.table(read.csv(infile_data, sep = ' '))
  data = data[data$Z < 0.5,]
  dists = comoving_distances_at_z(data$Z, cosmo$Om0, cosmo$OmK, cosmo$OmL, cosmo$H0)


  rho_mean = density_function(random_z, length(random_z)/400, g09_frac_area, cosmo)
  linking_lengths = rho_mean(data$Z)^(-1./3)
  max_on_sky_radius = coshaloMvirToRvir(Mmax, z = data$Z, Rho="crit", cosmo$H0)
  max_los_distances = coshaloMvirToSigma(Mmax, z= data$Z, Rho="crit", cosmo$H0)

  too_big = linking_lengths * b0 > max_on_sky_radius
  too_long = linking_lengths * b0 * r0 > max_los_distances

  linking_lengths[too_big] = max_on_sky_radius[too_big]/b0
  linking_lengths[too_long] = max_los_distances[too_long]/(b0 * r0)

  b0 = 0.06
  r0 = 18.
  rm(randoms)

  df = data.frame(ra = data$RA, dec = data$DEC, co_dists = dists, ll = linking_lengths)
  write.csv(df, "test_arrays.csv")

  #dumb_links = fof_links(df$ra, df$dec, df$co_dists, df$ll, b0, r0)
  groups = .find_groups(df$ra, df$dec, df$co_dists, df$ll, b0, r0)

})
