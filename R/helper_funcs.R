library(Highlander)
library(psych)

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
calculate_s_total <- function(measured_ids, group_ids, min_group_size=2) {
  return(calculate_s_score(as.integer(measured_ids), as.integer(group_ids), as.integer(min_group_size)))
}

#' Tuning the group finder to a mock catalogue to find the optimum b0 and r0.
#' @description
#' Uses Highlander(CMA+MCMC fitting) to find the optimum parameters for b0 and r0 given a mock catalogue
#' or list of mock catalgoues.
#'
#' @details
#' This function will average over all lightcones that are passed to it making it easy to pass
#' multiple independent lightcones which some users may wish to do to account for cosmic variance.
#' At the same time, the user can pass a single lightcone if that is all they wish. This function
#' is already parralized in multiple different ways so the user should not try and parralize it themselves
#' as this can actually slow the process down significantly. It should scale easily with number of cores.
#'
#' This is a helper function and it isn't necessary to use for tuning. Users can define the optimum
#' parameters in any way they wish. This is just a good approach that should meet most use cases.
#'
#' @param list_of_catalogs A list of RedshiftCatalog objects with mock_group_ids already set.
#' @param minimum_group_size The minim size of groups that should be accounted for when tuning.
#' @param b0_estimate The initial guess of the b0 linking length
#' @param r0_estimate The initial fiess of the r0 linking length
#' @param b0_bounds The lower and upper limits to search for the b0 linking length (e.g. c(0.03, 0.1)).
#' @param r0_bounds The lower and upper limits to search for the r0 linking length (e.g. c(5, 50)).
#' @param n_iterations The Number of iterations per CMA and Laplace's Demon respectively. Set to (100, 100)
#' @param n_final_mcmc The Number of iterations to run for the final MCMC run. Set to 2500
#' @returns A list of output from `Highlander` containg at least the best parameters of all iterations. See `Highlander` for more details.
#' @export
tune_group_finder <- function(
    list_of_catalogs,
    minimum_group_size,
    b0_estimate,
    r0_estimate,
    b0_bounds,
    r0_bounds,
    n_iterations = c(100, 100),
    nfinal_mcmc = 2500) {

  optimal_function <- function(par, redshift_catalogues) {
    b0 <- par[1]
    r0 <- par[2]
    message(cat(par))
    s_totals <- list()

    for (i in seq_along(redshift_catalogues)) {
      redshift_catalogue <- redshift_catalogues[[i]]

      # Skip if not a RedshiftCatalog
      if (!inherits(redshift_catalogue, "RedshiftCatalog")) {
        next
      }

      redshift_catalogue$run_fof(b0, r0)
      s_total <- redshift_catalogue$compare_to_mock(min_group_size = minimum_group_size)

      s_totals[[length(s_totals) + 1]] <- s_total

    }
    FoM <- psych::harmonic.mean(unlist(s_totals))
    if (is.nan(FoM)) {
      FoM <- 0
    }
    message('LP: ', FoM)
    message('-----------------')
    return(FoM)
  }

  opt_gama <- Highlander::Highlander(
    c(b0_estimate, r0_estimate),
    Data = list_of_catalogs,
    likefunc = optimal_function,
    likefunctype = "CMA",
    optim_iters = 2,
    liketype = "max",
    Niters = n_iterations,
    NfinalMCMC = nfinal_mcmc,
    lower = c(b0_bounds[1], r0_bounds[1]),
    upper = c(b0_bounds[2], r0_bounds[2]),
    parm.names = c("b0", "r0")
  )

  return(opt_gama)

}
