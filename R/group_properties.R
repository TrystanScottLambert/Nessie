
#' Generating the group catalog
#' @param ra The right ascension of the entire redshift catalog.
#' @param dec The declination of the entire redshift catalog.
#' @param redshift The redshifts of the entire redshift catalog.
#' @param magnidues The magnitudes for the entire redshift catalog.
#' @param velocity_erros The velocity errors for the entire redshift catalog.
#' @param group_ids The IDs of the groups for every galaxy in the redshift catalog.
#' @param cosmology a Flat LCDM cosmology object. This can be generated with `create_flat_cosmology`
#' @return Named list with group centers, radii, and multiplicity
generate_group_catalog <- function(ra, dec, redshift, magnitudes, velocity_errors, group_ids, cosmology) {
  create_group_catalog(
    ra,
    dec,
    redshift,
    magnitudes,
    velocity_errors,
    group_ids,
    cosmology$omega_m,
    cosmology$omega_k,
    cosmology$omega_l,
    cosmology$h0)
}
