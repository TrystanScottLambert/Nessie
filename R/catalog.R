library(R6)

RedshiftCatalog <- R6::R6Class("RedshiftCatalog",
  public = list(
    ra_array = NULL,
    dec_array = NULL,
    redshift_array = NULL,
    density_function = NULL,
    cosmology = NULL,
    completeness = NULL,
    group_ids = NULL,
    current_r0 = NULL,
    current_b0 = NULL,

    initialize = function(ra_array, dec_array, redshift_array, density_function, cosmology, completeness=NULL) {
      self$ra_array <- ra_array
      self$dec_array <- dec_array
      self$redshift_array <- redshift_array
      self$density_function <- density_function
      self$cosmology <- cosmology
      if (is.null(completeness)) {
        self$completeness <- rep(1, length(self$ra_array))
      } else {
        self$completeness <- completeness
      }
    },

    fof = function(b0, r0, max_stellar_mass = 1e15) {
      co_dists <- comoving_distances_at_z(self$redshift_array, self$cosmology$Om0, self$cosmology$OmK, self$cosmology$OmL, self$cosmology$H0)
      # Calculating the plane-of-sky linking lengths
      linking_lengths <- density_function(self$redshift_array)^(-1./3) * (self$completeness)^(-1./3)
      gal_rad <- b0 * linking_lengths
      max_on_sky_radius <- celestial::coshaloMvirToRvir(max_stellar_mass, z = self$redshift_array, Rho = "crit", self$cosmology$H0)
      too_wide <-  gal_rad > max_on_sky_radius
      gal_rad[too_wide] <- max_on_sky_radius[too_wide]
      linking_lengths_pos <- gal_rad/(self$cosmology$h * co_dists)

      # Calculating the line-of-sight linking lengths
      R <- r0 * (1 + self$redshift_array)/(sqrt(self$cosmology$Om0 * (1 + self$redshift_array)^3 + self$cosmology$OmL))
      linking_lengths_los <- (gal_rad * R)/self$cosmology$h
      max_los_distances <- celestial::coshaloMvirToSigma(max_stellar_mass, z = self$redshift_array, Rho = "crit", self$cosmology$H0) * (1 + self$redshift_array) / celestial::cosgrowH(H0 = self$cosmology$H0, z = self$redshift_array)
      too_far <- linking_lengths_los > max_los_distances
      linking_lengths_los[too_far] <- max_los_distances[too_far]
      groups <- .find_groups(self$ra_array, self$dec_array, co_dists, linking_lengths_pos, linking_lengths_los)
      return(groups)
    },

    set_group_ids = function(b0, r0, max_stellar_mass = 1e15) {
      group_links <- self$fof(b0, r0, max_stellar_mass)
      all_ids <- seq(length(self$ra_array)) # positions of all galaxies
      singleton_galaxies <- setdiff(all_ids, group_links$galaxy_id)
      singleton_marker_id <- rep(-1, length(singleton_galaxies))
      all_galaxies <- rbind(group_links, data.frame(galaxy_id = singleton_galaxies, group_id = singleton_marker_id))
      self$group_ids <- as.integer(all_galaxies[order(all_galaxies$galaxy_id), ]$group_id)
    },

    calculate_group_table = function(absolute_magnitudes, velocity_errors) {
      as.data.frame(create_group_catalog(
        self$ra_array, self$dec_array, self$redshift_array, absolute_magnitudes, velocity_errors,
        self$group_ids, self$cosmology$Om0, self$cosmology$OmK, self$cosmology$OmL,
        self$cosmology$H0))
    }
  )
)
