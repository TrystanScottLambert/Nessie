#' Running density function estimation
#'
#' @description
#' `create_density_function` returns the rho(z) function (Eq ?? in Lambert+2025).
#' This is a required argument for the FoF algorithm
#'
#' @details
#' The n(z) is approximated from a given array of redshift values. These values can be generated
#' in any way. In the simpliest case this can just be the redshifts of the redshift survey.
#' More approriate methods should be used though. The Probability Distribution function is then
#' estimated and integrated over several bins. The total running density is simply the total counts
#' in a bin multipled by the PDF and divided by the running volume.
#'
#' @param redshifts An array of redshift values.
#' @param total_counts The total counts of the redshift redshift survey. Doesn't make sense.
#' @param survey_fractional_area The fraction of the survey over a full 4pi steradian area.
#' @param cosmology A FlatCosmology object (created using the FlatCosmology class).
#' @param binwidth The binwidth (in comoving distance) which will be used for the integration
#' @param N The total number of data points to interpolate over.
#'
#' @return running_density A function which is the rho(z).
#' @export
create_density_function <- function(redshifts, total_counts, survey_fractional_area, cosmology, binwidth=40, N = 1e4) {

  comoving_distances <- cosmology$comoving_distance(redshifts)
  max_comoving <- max(comoving_distances) + 2*binwidth # allowing a 2 binwidth grace
  #TODO: Remove this shit when before release
  max_comoving <- 2000
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
  running_density_z <- approxfun(cosmology$z_at_comoving_dist(kde$x), (total_counts * running_integral)/running_volume, rule = 2)
  return(running_density_z)
}


#' Comparing measured groups to mock catalogue groups
#' @description
#' Compares the current group_ids to a mock known grouping ids.
#' @param mock_group_ids group ids from the mock catalog. Have to be integers.
#' @param singleton_ids the integer id value that is assigned to all singleton galaxies.
#' @returns a list containing all the values from equations 9 - 15 in Robotham+2011
#' @export
calculate_mock_comparison_metrics <- function(measured_ids, group_ids, singleton_id) {
  return(calculate_cost_metrics(as.integer(measured_ids), as.integer(group_ids), as.integer(singleton_id)))
}
