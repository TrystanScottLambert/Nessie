
#' Create a Flat LCDM cosmology parameter object.
#' @description
#' `create_flat_cosmology` will return a named list with the cosmological parameters
#' with OmL = 1 - Om0.
#' @param h Little h. Often taken to be 0.7. H0 = 100 * h.
#' @param omega_m Density of matter in the universe. Usually taken to be 0.3
#' @return parameters A named list with the cosmological parameters.
#' @examples
#' cosmo <- create_cosmology(H0 = 70, Om0 = 0.3)
#' cosmo$H0  # Access via $
create_flat_cosmology <- function(h, omega_m) {
  parameters <- list(
    h = h,
    H0 = 100 * h,
    Om0 = omega_m,
    OmL = 1 - omega_m,
    Omk = 0.,
    OmR = 0.
  )
  class(parameters) <- "cosmology"
  return(parameters)
}


#' Running density function estimation
#'
#' @description
#' `density_function` returns the rho(z) function (Eq ?? in Lambert+2025).
#' This is a required argument for the FoF algoritm
#' 
#' @details
#' The n(z) is approximated from a given array of redshift values. These values can be generated
#' in anyway. In the simpliest case this can just be the redshifts of the redshift survey. 
#' More approriate methods should be used though. The Probability Distribution function is then
#' estimated and integrated over several bins. The total running density is simply the total counts
#' in a bin multipled by the PDF and divided by the running volume.
#'
#' @param redshifts An array of redshift values.
#' @param total_counts The total counts of the redshift redshift survey. Doesn't make sense.
#' @param survey_fractional_area The fraction of the survey over a full 4pi steradian area.
#' @param binwidth The binwidth (in comoving distance) which will be used for the integration
#' @param N The total number of data points to interpolate over.
#' 
#' @return running_density A function which is the rho(z). 
density_function <- function(redshifts, total_counts, survey_fractional_area, cosmology, binwidth=40, N = 1e4) {

    comoving_distances <- comoving_distances_at_z(redshifts, cosmology$Om0, cosmology$Omk, cosmology$OmL, cosmology$H0)
    max_comoving <- max(comoving_distances) + 2*binwidth # allowing a 2 bindwidth grace
    kde <- density(comoving_distances, bw = binwidth / sqrt(12), from = 0, to = max_comoving, n = N, kern = "rect")
    kde_func <- approxfun(kde$x, kde$y, rule = 2)

    comoving_bins <- seq(0, max_comoving, len = N)
    running_integral <- {}
    for (colim in comoving_bins) {
    running_integral <- c(running_integral, integrate(kde_func, colim - binwidth / 2, colim + binwidth / 2)$value)
    }

    upper_volumes <- (4 / 3) * pi * (comoving_bins + binwidth / 2)^3
    lower_volumes <- (4 / 3) * pi * (comoving_bins - binwidth / 2)^3
    running_volume <- survey_fractional_area * (upper_volumes - lower_volumes)

    running_density_z <- approxfun(z_at_comoving_distances(kde$x, cosmology), (total_counts * running_integral)/running_volume, rule=2)
    return(running_density_z)
}
