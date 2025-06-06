library(R6)

FlatCosmology <- R6::R6Class("FlatCosmology",
   public = list(
     h = NULL,
     hubble_constant = NULL,
     omega_m = NULL,
     omega_lambda = NULL,
     omega_k = NULL,
     omega_radiation = NULL,

     initialize = function(h, omega_matter) {
       self$h <- h
       self$hubble_constant <- 100 * self$h
       self$omega_m <- omega_matter
       self$omega_lambda <- 1 - self$omega_m
       self$omega_k <- 0.
       self$omega_radiation <- 0.
     },

     comoving_distance = function(redshift) {
       return(comoving_distances_at_z(redshift, self$omega_m, self$omega_k, self$omega_lambda, self$hubble_constant))
     },

     z_at_comoving_dist = function(co_dist) {
        return(z_at_comoving_distances(co_dist, self$omega_m, self$omega_k, self$omega_lambda, self$hubble_constant))
     },

     virial_radius = function(mass_solar, redshift) {
       celestial::coshaloMvirToRvir(mass_solar, z = redshift, Rho = "crit", H0 = self$hubble_constant, OmegaM = self$omega_m)
     },

     velocity_dispersion = function(mass_solar, redshift) {
       celestial::coshaloMvirToSigma(mass_solar, z = redshift_array, Rho = "crit", H0 = self$hubble_constant, OmegaM = self$omega_m)
     },

     h0_grow = function(redshift) {
       celestial::cosgrowH(H0 = self$hubble_constant, z = redshift, OmegaM = self$omega_m)
     }
  )
)


