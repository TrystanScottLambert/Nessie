% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_funcs.R
\name{.find_groups}
\alias{.find_groups}
\title{Finds groups in the given redshift arrays.}
\usage{
.find_groups(
  ra_array,
  dec_array,
  comoving_distances,
  linking_lengths_pos,
  linking_lengths_los
)
}
\arguments{
\item{ra_array}{An array-like object of right ascension in decimal degrees.}

\item{dec_array}{An array-like object of declination values in decimal degrees.}

\item{comoving_distances}{An array-like object of comoving_distances in Mpc.}

\item{linking_lengths_pos}{An array of the plane-of-sky linking lengths.}

\item{linking_lengths_los}{An array of the line-of-sight linking lengths.}
}
\value{
Data Frame of the galaxy ids and which groups they are in.
}
\description{
Identifies the groups based on the given data and linking length information.
}
\details{
This is a helper function which first identifies the pairs of galaxies (all the friends) in the
given data and then constructs a graph and identifies the connected components from these
pairs. This is done using two functions: \code{fof_links} and \code{.group_graph}. This function shouldn't
be used by the user, who are better served using the \code{fof} function directly.
}
