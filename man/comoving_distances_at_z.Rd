% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extendr-wrappers.R
\name{comoving_distances_at_z}
\alias{comoving_distances_at_z}
\title{Calculates multiple comoving distances for multiple redshifts.}
\usage{
comoving_distances_at_z(redshift_array, omega_m, omega_k, omega_l, h0)
}
\arguments{
\item{redshift_array}{an array of multiple redshift values.}

\item{omega_m}{Mass density (often 0.3 in LCDM).}

\item{omega_k}{Effective mass density of relativistic particles (often 0. in LCDM).}

\item{omega_l}{Effective mass density of dark energy (often 0.7 in LCDM).}

\item{h0}{H0 = 100 * h.}
}
\value{
multiple comoving distance in Mpc.
}
\description{
Calculates multiple comoving distances for multiple redshifts.
}
