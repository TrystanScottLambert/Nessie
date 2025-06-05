#' Create a Flat LCDM cosmology parameter object.
#' @description
#' `create_flat_cosmology` will return a named list with the cosmological parameters
#' with OmL = 1 - Om0.
#' @param h Little h. Often taken to be 0.7. H0 = 100 * h.
#' @param omega_m Density of matter in the universe. Usually taken to be 0.3
#' @returns parameters A named list with the cosmological parameters.
#' @examples
#' cosmo <- create_cosmology(H0 = 70, Om0 = 0.3)
#' cosmo$H0  # Access via $
create_flat_cosmology <- function(h, omega_m) {
  parameters <- list(
    h = h,
    H0 = 100 * h,
    Om0 = omega_m,
    OmL = 1 - omega_m,
    OmK = 0.,
    OmR = 0.
  )
  class(parameters) <- "cosmology"
  return(parameters)
}

#' Running density function estimation
#'
#' @description
#' `create_density_function` returns the rho(z) function (Eq ?? in Lambert+2025).
#' This is a required argument for the FoF algoritm
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
#' @param binwidth The binwidth (in comoving distance) which will be used for the integration
#' @param N The total number of data points to interpolate over.
#'
#' @return running_density A function which is the rho(z).
create_density_function <- function(redshifts, total_counts, survey_fractional_area, cosmology, binwidth=40, N = 1e4) {

    comoving_distances <- comoving_distances_at_z(redshifts, cosmology$Om0, cosmology$OmK, cosmology$OmL, cosmology$H0)
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
    running_density_z <- approxfun(z_at_comoving_distances(kde$x, cosmology$Om0, cosmology$OmK, cosmology$OmL, cosmology$H0), (total_counts * running_integral)/running_volume, rule=2)
    return(running_density_z)
}

#' Finding groups from a graph structure.
#'
#' @description
#' Creates a graph structure from the galaxy pairs and then finds the groups within those pairs.
#'
#' @details
#' The pairs of galaxies can be though of as edges of a graph. An undirected graph data structure
#' can be built from the edges. All the connected components of the graph would then be the groups.
#' These are retured as such.
#'
#' `.group_graph` should not be used by the user. This is a helper function for the fof function.
#'
#' @param links A named list with i and j. These links are made using the `fof_links` function.
#' @return A dataframe with galaxy ids in one column and group id in another.
#'
.group_graph <-  function(links){
  group = igraph::components(igraph::graph_from_data_frame(links, directed=FALSE))
  group = data.frame(galaxy_id = as.integer(names(group$membership)), group_id=group$membership)
  return(group)
}

#' Finds groups in the given redshift arrays.
#'
#' @description
#' Identifies the groups based on the given data and linking length information.
#'
#' @details
#' This is a helper function which first identifies the pairs of galaxies (all the friends) in the
#' given data and then constructs a graph and identifies the connected compoenents from these
#' pairs. This is done using two functions: `fof_links` and `.group_graph`. This function shouldn't
#' be used by the user, who are better served using the `fof` function directly.
#'
#' @param ra_array An array-like object of right ascension in decimal degrees.
#' @param dec_array An array-like object of declination values in decimal degrees.
#' @param comoving_distances An array-like object of comoving_distances in Mpc.
#' @param linking_lengths An array-like object of linking lengths for every galaxy. These are before
#' accounting for the constants, b0, and r0.
#' @param b0 Plane-of-sky constant which will be scaled by the given linking lengths.
#' @param r0 Line-of_sight constant value which will be scaled by the given linking lengths and b0.
#'
#' @return Data Frame of the galaxy ids and which groups they are in.
#'
.find_groups <- function(ra_array, dec_array, comoving_distances, linking_lengths_pos, linking_lengths_los) {
  links <- fof_links_fast(ra_array, dec_array, comoving_distances,  linking_lengths_pos, linking_lengths_los)
  group <- .group_graph(links)
  return(group)
}

.find_groups_aaron <- function(ra_array, dec_array, comoving_distances, linking_lengths_pos, linking_lengths_los) {
  links <- fof_links_aaron(ra_array, dec_array, comoving_distances,  linking_lengths_pos, linking_lengths_los)
  group <- .group_graph(links)
  return(group)
}

#' Run the friends-of-friends algorithm
#'
#' @description
#' Performs the friends-of-friends algorithm on the given data.
#'
#' @details
#' This is the core function of the group finder and actually runs the friends-of-friends algorithm.
#' Before running this function it is important to have the density function already calculated.
#' This can be done by passing an array of n(z) values to the `density_function`. An array of completeness
#' can also be passed if needed or else 100% completeness is assumed.
#'
#' @param b0 Plane-of-sky constant which will be scaled by the given linking lengths.
#' @param r0 Line-of_sight constant value which will be scaled by the given linking lengths and b0.
#' @param ras An array-like object of right ascension in decimal degrees.
#' @param decs An array-like object of declination values in decimal degrees.
#' @param redshifts An array-like object of comoving_distances in Mpc.
#' @param density_function An array-like object of linking lengths for every galaxy. These are before
#' @param cosmology A cosmology object. (can be created with `create_flat_cosmology`).
#' @param completeness An array of equal size to the data with values between 0 and 1. Representing how
#' complete the survey is around that particular galaxy. This will scale that galaxy's linking lengths
#' accordingly.
#'
#' @return Data Frame of the galaxy ids and which groups they are in.
#'
fof <- function(b0, r0, ras, decs, redshifts, density_function, cosmology, completeness = rep(1, length(ras)), Mmax=1e15) {
  co_dists <- comoving_distances_at_z(redshifts, cosmology$Om0, cosmology$OmK, cosmology$OmL, cosmology$H0)

  # Calculating the plane-of-sky linking lengths
  linking_lengths = density_function(redshifts)^(-1./3) * (completeness)^(-1./3)
  gal_rad = b0 * linking_lengths
  max_on_sky_radius = celestial::coshaloMvirToRvir(Mmax, z = redshifts, Rho="crit", cosmology$H0)
  too_wide =  gal_rad > max_on_sky_radius
  gal_rad[too_wide] = max_on_sky_radius[too_wide]
  linking_lengths_pos = gal_rad/(cosmology$h * co_dists)

  # Calculating the line-of-sight linking lengths
  R = r0 * (1 + redshifts)/(sqrt(cosmology$Om0 * (1+redshifts)^3 + cosmology$OmL))
  linking_lengths_los = (gal_rad * R)/cosmology$h
  max_los_distances = celestial::coshaloMvirToSigma(Mmax, z= redshifts, Rho="crit", cosmology$H0) * (1+redshifts) / celestial::cosgrowH(H0=cosmology$H0, z=redshifts)
  too_far = linking_lengths_los > max_los_distances
  linking_lengths_los[too_far] = max_los_distances[too_far]

  return(.find_groups(ras, decs, co_dists, linking_lengths_pos, linking_lengths_los))
}

.fof_aaron <- function(b0, r0, ras, decs, redshifts, density_function, cosmology, completeness = rep(1, length(ras)), Mmax=1e15) {
  co_dists <- comoving_distances_at_z(redshifts, cosmology$Om0, cosmology$OmK, cosmology$OmL, cosmology$H0)

  # Calculating the plane-of-sky linking lengths
  linking_lengths = density_function(redshifts)^(-1./3) * (completeness)^(-1./3)
  gal_rad = b0 * linking_lengths
  max_on_sky_radius = celestial::coshaloMvirToRvir(Mmax, z = redshifts, Rho="crit", cosmology$H0)
  too_wide =  gal_rad > max_on_sky_radius
  gal_rad[too_wide] = max_on_sky_radius[too_wide]
  linking_lengths_pos = gal_rad/(cosmology$h * co_dists)

  # Calculating the line-of-sight linking lengths
  R = r0 * (1 + redshifts)/(sqrt(cosmology$Om0 * (1+redshifts)^3 + cosmology$OmL))
  linking_lengths_los = (gal_rad * R)/cosmology$h
  max_los_distances = celestial::coshaloMvirToSigma(Mmax, z= redshifts, Rho="crit", cosmology$H0) * (1+redshifts) / celestial::cosgrowH(H0=cosmology$H0, z=redshifts)
  too_far = linking_lengths_los > max_los_distances
  linking_lengths_los[too_far] = max_los_distances[too_far]

  return(.find_groups_aaron(ras, decs, co_dists, linking_lengths_pos, linking_lengths_los))
}



#' Get the group ids for the given redshift survey.
#'
#' @description
#' Performs the FoF algorithm on the given data and returns the group_id column for that redshift survey.
#'
#' @details
#' This function is a helper function which wraps around the `fof` function but instead of returning
#' the galaxy_id and associated group, this function will return a single array of group_ids which is
#' equal to the length of the arrays. The ids here do not differentiate between binaries and groups.
#' Galaxies which were not found in a group are given the -1 group id.
#'
#' This is built to work well with the `generate_group_catalog` function and will put the group ids
#' in the order that is required for that function.
#'
#' @param b0 Plane-of-sky constant which will be scaled by the given linking lengths.
#' @param r0 Line-of_sight constant value which will be scaled by the given linking lengths and b0.
#' @param ras An array-like object of right ascension in decimal degrees.
#' @param decs An array-like object of declination values in decimal degrees.
#' @param redshifts An array-like object of comoving_distances in Mpc.
#' @param density_function An array-like object of linking lengths for every galaxy. These are before
#' @param cosmology A cosmology object. (can be created with `create_flat_cosmology`).
#' @param completeness An array of equal size to the data with values between 0 and 1. Representing how
#' complete the survey is around that particular galaxy. This will scale that galaxy's linking lengths
#' accordingly.
#'
#' @return A single array of group ids equal to the length of the arrays that were given.
#'
get_group_ids <- function(b0, r0, ras, decs, redshifts, density_function, cosmology, completeness = rep(1, length(ras))) {
  group_links <- fof(b0, r0, ras, decs, redshifts, density_function, cosmology, completeness)
  all_ids <- seq(length(ras)) # positions of all galaxies
  singleton_galaxies <- setdiff(all_ids, group_links$galaxy_id)
  singleton_marker_id <- rep(-1, length(singleton_galaxies))
  all_galaxies = rbind(group_links, data.frame(galaxy_id = singleton_galaxies, group_id = singleton_marker_id))
  return(all_galaxies[order(all_galaxies$galaxy_id), ]$group_id)
}



