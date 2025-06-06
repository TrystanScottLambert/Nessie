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
  group = igraph::components(igraph::graph_from_data_frame(links, directed = FALSE))
  group = data.frame(galaxy_id = as.integer(names(group$membership)), group_id = group$membership)
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

.find_groups_classic <- function(ra_array, dec_array, comoving_distances, linking_lengths_pos, linking_lengths_los) {
  links <- fof_links_aaron(ra_array, dec_array, comoving_distances,  linking_lengths_pos, linking_lengths_los)
  group <- .group_graph(links)
  return(group)
}

