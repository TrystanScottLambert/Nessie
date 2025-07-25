% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extendr-wrappers.R
\name{fof_links_aaron}
\alias{fof_links_aaron}
\title{finding the links between all galaxies in a brute force way.}
\usage{
fof_links_aaron(
  ra_array,
  dec_array,
  comoving_distances,
  linking_lengths_pos,
  linking_lengths_los
)
}
\arguments{
\item{comoving_distances}{Array of comoving distances in Mpc.}

\item{ra}{Array of right ascension values.}

\item{dec}{Array of declination values.}

\item{linking_lengths}{An array of individual scaled linking lengths for each galaxy (ignoring r0 and b0).}

\item{b0}{The plane-of-sky constant to be scaled.}

\item{r0}{The line-of-sight constant to be scaled.}
}
\value{
A dataframe-like object of tuples which represent the link between galaxies (i, j) if they exist.
}
\description{
\code{fof_links_aaron} will determine all connections between galaxies in a survey and return the pairs.
}
