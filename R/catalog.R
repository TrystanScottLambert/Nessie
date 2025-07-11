library(R6)

#' RedshiftCatalog class
#'
#' @description
#' A catalog object for redshift-space galaxy data that supports group finding
#' using a Friends-of-Friends (FoF) algorithm with redshift-dependent linking lengths.
#'
#' @field ra_array A numeric vector of Right Ascension values (in degrees).
#' @field dec_array A numeric vector of Declination values (in degrees).
#' @field redshift_array A numeric vector of redshift values.
#' @field density_function A function that takes redshift values and returns the local density.
#' @field cosmology A FlatCosmology object (created using the FlatCosmology class).
#' @field completeness A numeric vector of completeness values per galaxy (default is 1 for all).
#' @field group_ids An integer vector storing the resulting group ID for each galaxy (after `run_fof()` is called).
#' @field current_r0 The most recent `r0` used in `run_fof()`
#' @field current_b0 The most recent `b0` used in `run_fof()`
#' @field mock_group_ids The IDs of 'true' groupings from a mock catalogue in the same format at group_ids. Singleton groups need an id of -1. (only used when tuning).
#'
#' @export
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
    mock_group_ids = NULL,

    #' @description
    #' Create a new RedshiftCatalog object.
    #' @param ra_array A numeric vector of right ascension values (degrees).
    #' @param dec_array A numeric vector of declination values (degrees).
    #' @param redshift_array A numeric vector of redshifts.
    #' @param density_function A function that takes redshifts and returns galaxy density.
    #' @param cosmology A FlatCosmology object (created using the FlatCosmology class).
    #' @param completeness An optional numeric vector of completeness weights (default: 1 for all galaxies).
    initialize = function(
      ra_array, dec_array, redshift_array, density_function, cosmology, completeness=NULL
      ) {
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


    #' @description
    #' Generate FoF links between galaxies based on spatial and redshift linking lengths. There is
    #' very little reason to have to run this yourself. In most cases it is more appropriate to run
    #' the `run_fof` method.
    #' @param b0 The plane-of-sky linking length constant.
    #' @param r0 The line-of-sight linking length constant.
    #' @param max_stellar_mass The maximum stellar mass to cap linking distances (default: 1e15 solar masses).
    #' @return A data.frame with columns `galaxy_id` and `group_id` indicating linked galaxies.
    get_raw_groups = function(b0, r0, max_stellar_mass = 1e15, algorithm = 'fast') {
      co_dists <- self$cosmology$comoving_distance(self$redshift_array)
      # Calculating the plane-of-sky linking lengths
      linking_lengths <- self$density_function(self$redshift_array)^(-1./3) * (self$completeness)^(-1./3)
      gal_rad <- b0 * linking_lengths
      max_on_sky_radius <- self$cosmology$virial_radius(max_stellar_mass, self$redshift_array)
      too_wide <-  gal_rad > max_on_sky_radius
      gal_rad[too_wide] <- max_on_sky_radius[too_wide]
      linking_lengths_pos <- gal_rad/(self$cosmology$h * co_dists)

      # Calculating the line-of-sight linking lengths
      R <- r0 * (1 + self$redshift_array) /
        (sqrt(self$cosmology$omega_m * (1 + self$redshift_array)^3 + self$cosmology$omega_lambda))
      linking_lengths_los <- (gal_rad * R)/self$cosmology$h
      max_los_distances <- self$cosmology$velocity_dispersion(max_stellar_mass, self$redshift_array) * (1 + self$redshift_array) / self$cosmology$h0_grow(self$redshift_array)
      too_far <- linking_lengths_los > max_los_distances
      linking_lengths_los[too_far] <- max_los_distances[too_far]
      if (algorithm == 'fast') {
        groups <- .find_groups(self$ra_array, self$dec_array, co_dists, linking_lengths_pos, linking_lengths_los)
      } else {
        groups <- .find_groups_classic(self$ra_array, self$dec_array, co_dists, linking_lengths_pos, linking_lengths_los)
      }
      return(groups)
    },

    #' @description
    #' Run the full Friends-of-Friends (FoF) algorithm and assign group IDs to all galaxies.
    #' Singleton galaxies (unlinked) are given group ID -1.
    #' @param b0 The plane-of-sky linking length scale factor.
    #' @param r0 The line-of-sight linking length scale factor.
    #' @param max_stellar_mass The maximum stellar mass to cap linking distances (default: 1e15).
    run_fof = function(b0, r0, max_stellar_mass = 1e15) {
      group_links <- self$get_raw_groups(b0, r0, max_stellar_mass)
      all_ids <- seq(length(self$ra_array)) # positions of all galaxies
      singleton_galaxies <- setdiff(all_ids, group_links$galaxy_id)
      singleton_marker_id <- rep(-1, length(singleton_galaxies))
      all_galaxies <- rbind(
        group_links, data.frame(galaxy_id = singleton_galaxies, group_id = singleton_marker_id))
      self$group_ids <- as.integer(all_galaxies[order(all_galaxies$galaxy_id), ]$group_id)
      self$current_r0 <- r0
      self$current_b0 <- b0
    },

    #' @description
    #' Generate a summary data.frame of group properties based on assigned group IDs. Must have run the group finder.
    #' @param absolute_magnitudes A numeric vector of absolute magnitudes per galaxy.
    #' @param velocity_errors A numeric vector of redshift or velocity errors.
    #' @return A data.frame summarizing group-level statistics.
    calculate_group_table = function(absolute_magnitudes, velocity_errors) {
      if (is.null(self$group_ids)) {
        stop("No group ids found. Be sure to run the `run_fof` method before calling `calculate_group_table`")
      }
      as.data.frame(create_group_catalog(
        self$ra_array, self$dec_array, self$redshift_array, absolute_magnitudes, velocity_errors,
        self$group_ids, self$cosmology$omega_m, self$cosmology$omega_k, self$cosmology$omega_lambda,
        self$cosmology$hubble_constant))
    },

    #' @description
    #' Compares the current group_ids to a mock known grouping ids. Must have run the group finder and set both the mock_group_ids and singleton_id
    #' @details
        #' This is a wrapper around the `calculate_mock_comparison_metrics` function which is available in Nessie.
        #' But this is designed to work specifically within the RedshiftCatalog class.
        #' `run_fof` should be run first. The `mock_group_ids` and `singleton_id` should be set (red_cat$mock_groups_ids = mock_group_ids)
    #' @param min_group_size The minimum size of a group that should be included in the cost metric.
    #' @returns a list containing all the values from equations 9 - 15 in Robotham+2011
    compare_to_mock = function(min_group_size = 2) {
      if (is.null(self$group_ids)) {
        stop("No group ids found. Be sure to run the `run_fof` method before calling `calculate_group_table`")
      }

      if (is.null(self$mock_group_ids)) {
        stop("No mock group ids found. Be sure to set the comparison mock group id with `$mock_group_ids = ` ")
      }

      return(calculate_s_total(self$group_ids, self$mock_group_ids, min_group_size))
    }
  )
)
