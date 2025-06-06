#' Finding groups from a graph structure.
#'
#' @description
#' Creates a graph structure from the galaxy pairs and then finds the groups within those pairs.
#'
#' @details
#' The pairs of galaxies can be though of as edges of a graph. An undirected graph data structure
#' can be built from the edges. All the connected components of the graph would then be the groups.
#' These are returned as such.
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
#' given data and then constructs a graph and identifies the connected components from these
#' pairs. This is done using two functions: `fof_links` and `.group_graph`. This function shouldn't
#' be used by the user, who are better served using the `fof` function directly.
#'
#' @param ra_array An array-like object of right ascension in decimal degrees.
#' @param dec_array An array-like object of declination values in decimal degrees.
#' @param comoving_distances An array-like object of comoving_distances in Mpc.
#' @param linking_lengths_pos An array of the plane-of-sky linking lengths.
#' @param linking_lengths_los An array of the line-of-sight linking lengths.
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

