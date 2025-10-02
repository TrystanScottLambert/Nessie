library(R6)

#' FlatCosmology class
#'
#' @description
#' An R6 class for flat ΛCDM cosmology models, providing methods for computing
#' comoving distances, virial radii, velocity dispersions, and Hubble evolution.
#'
#' @field h Dimensionless Hubble parameter (i.e., H0 / 100).
#' @field hubble_constant Hubble constant in km/s/Mpc.
#' @field omega_m Matter density parameter (Ωₘ).
#' @field omega_lambda Cosmological constant / dark energy density (Ω_Λ), derived assuming flatness.
#' @field omega_k Curvature density parameter (Ωₖ), fixed to 0 in flat cosmology.
#' @field omega_radiation Radiation density parameter (Ωᵣ), fixed to 0 by default.
#'
#' @export
FlatCosmology <- R6::R6Class("FlatCosmology",
  public = list(
    h = NULL,
    hubble_constant = NULL,
    omega_m = NULL,
    omega_lambda = NULL,
    omega_k = NULL,
    omega_radiation = NULL,

    #' @description
    #' Create a new FlatCosmology object.
    #' @param h Dimensionless Hubble parameter (H0 / 100).
    #' @param omega_matter Matter density parameter Ωₘ.
    initialize = function(h, omega_matter) {
      self$h <- h
      self$hubble_constant <- 100 * self$h
      self$omega_m <- omega_matter
      self$omega_lambda <- 1 - self$omega_m
      self$omega_k <- 0.
      self$omega_radiation <- 0.
    },

    #' @description
    #' Compute comoving distances for a vector of redshifts.
    #' @param redshift A numeric vector of redshifts.
    #' @return A numeric vector of comoving distances (in Mpc).
    comoving_distance = function(redshift) {
      return(comoving_distances_at_z(redshift, self$omega_m, self$omega_k, self$omega_lambda, self$hubble_constant))
    },

    #' @description
    #' Compute the distance modulus for a vector of redshifts.
    #' @param redshift A numeric vector of redshifts.
    #' @return A numeric vector of distance moduli.
    dist_mod = function(redshift) {
      return(distance_modulus(redshift, self$omega_m, self$omega_k, self$omega_lambda, self$hubble_constant))
    },

    #' @description
    #' Compute redshift values corresponding to comoving distances.
    #' @param co_dist A numeric vector of comoving distances (in Mpc).
    #' @return A numeric vector of redshifts.
    z_at_comoving_dist = function(co_dist) {
      return(z_at_comoving_distances(co_dist, self$omega_m, self$omega_k, self$omega_lambda, self$hubble_constant))
    },

    #' @description
    #' Compute the virial radius of halos with given mass and redshift.
    #' @param mass_solar A numeric vector of halo masses in solar masses.
    #' @param redshift A numeric vector of redshifts.
    #' @return A numeric vector of virial radii (in kpc).
    virial_radius = function(mass_solar, redshift) {
      calculate_max_rvirs(mass_solar, redshift, self$omega_m, self$omega_k, self$omega_lambda, self$hubble_constant)
    },

    #' @description
    #' Compute the velocity dispersion of halos with given mass and redshift.
    #' @param mass_solar A numeric vector of halo masses in solar masses.
    #' @param redshift A numeric vector of redshifts.
    #' @return A numeric vector of velocity dispersions (in km/s).
    velocity_dispersion = function(mass_solar, redshift) {
      calculate_max_sigmas(mass_solar, redshift, self$omega_m, self$omega_k, self$omega_lambda, self$hubble_constant)
    },

    #' @description
    #' Compute the redshift-dependent Hubble parameter H(z) in a flat cosmology.
    #' @param redshift A numeric vector of redshifts.
    #' @return A numeric vector of H(z) values (in km/s/Mpc).
    h0_grow = function(redshift) {
      h_at_z(redshift, self$omega_m, self$omega_k, self$omega_lambda, self$hubble_constant)
    }

    #' @description
    #' Compute the redshift-dependent Hubble parameter H(z) in a flat cosmology.
    #' @param redshift A numeric vector of redshifts.
    #' @return A numeric vector of differential comoving volume in Mpc/sr.
    differential_covol = function(redshift) {
      diff_covol(redshift, self$omega_k, self$omega_lambda, self$hubble_constant)
    }

  )
)
