use extendr_api::prelude::*;
use fof::bijectivity::s_score;
use fof::completeness::{calculate_completeness, PositionCatalog};
use fof::group_properties::GroupedGalaxyCatalog;
use fof::link_finder::{ffl1, find_links};
use fof::stats::harmonic_mean;
use fof::Cosmology;
use randoms::generate_randoms;
use rayon::prelude::*;
use std::f64;

/// Calculate the hubble constant at different redshifts.
/// @param redshift_array an array of multiple redshift values.
/// @param omega_m Mass density (often 0.3 in LCDM).
/// @param omega_k Effective mass density of relativistic particles (often 0. in LCDM).
/// @param omega_l Effective mass density of dark energy (often 0.7 in LCDM).
/// @param h0 H0 = 100 * h.
/// @returns Multiple H(z) for different z.
/// @export
#[extendr]
fn h_at_z(redshift_array: Vec<f64>, omega_m: f64, omega_k: f64, omega_l: f64, h0: f64) -> Vec<f64> {
    let cosmo = Cosmology {
        omega_m,
        omega_k,
        omega_l,
        h0,
    };
    redshift_array
        .par_iter()
        .map(|z| cosmo.h_at_z(*z))
        .collect()
}

/// Generates random redshifts to model the n(z) without underlying LSS.
/// @param redshifts Redshifts from the survey.
/// @param mags Magnitudes from the survey.
/// @param z_lim The maximum z over which to confine all calculations.
/// @param maglim The magnitude limit of the survey.
/// @param n_clone The number of times more the randoms catalogue should be
/// @param iterations The number of iterations to iteratively iterate. (5-10 should be plenty)
/// @param omega_m Mass density (often 0.3 in LCDM).
/// @param omega_k Effective mass density of relativistic particles (often 0. in LCDM).
/// @param omega_l Effective mass density of dark energy (often 0.7 in LCDM).
/// @param h0 H0 = 100 * h.
#[extendr]
fn gen_randoms(
    redshifts: Vec<f64>,
    mags: Vec<f64>,
    z_lim: f64,
    maglim: f64,
    n_clone: i32,
    iterations: i32,
    omega_m: f64,
    omega_k: f64,
    omega_l: f64,
    h0: f64,
) -> Vec<f64> {
    let cosmo = randoms::cosmology::Cosmology {
        omega_m,
        omega_k,
        omega_l,
        h0,
    };
    generate_randoms(redshifts, mags, z_lim, maglim, n_clone, iterations, cosmo)
}

/// Calculates multiple comoving distances for multiple redshifts.
/// @param redshift_array an array of multiple redshift values.
/// @param omega_m Mass density (often 0.3 in LCDM).
/// @param omega_k Effective mass density of relativistic particles (often 0. in LCDM).
/// @param omega_l Effective mass density of dark energy (often 0.7 in LCDM).
/// @param h0 H0 = 100 * h.
/// @return multiple comoving distance in Mpc.
/// @export
#[extendr]
fn comoving_distances_at_z(
    redshift_array: Vec<f64>,
    omega_m: f64,
    omega_k: f64,
    omega_l: f64,
    h0: f64,
) -> Vec<f64> {
    let cosmo = Cosmology {
        omega_m,
        omega_k,
        omega_l,
        h0,
    };
    redshift_array
        .par_iter()
        .map(|z| cosmo.comoving_distance(*z))
        .collect()
}

/// Redshift at some given comoving distances in Mpc.
/// @param redshift_array an array of multiple redshift values.
/// @param omega_m Mass density (often 0.3 in LCDM).
/// @param omega_k Effective mass density of relativistic particles (often 0. in LCDM).
/// @param omega_l Effective mass density of dark energy (often 0.7 in LCDM).
/// @param h0 H0 = 100 * h.
/// @export
#[extendr]
fn z_at_comoving_distances(
    distances: Vec<f64>,
    omega_m: f64,
    omega_k: f64,
    omega_l: f64,
    h0: f64,
) -> Vec<f64> {
    let cosmo = Cosmology {
        omega_m,
        omega_k,
        omega_l,
        h0,
    };
    distances
        .par_iter()
        .map(|d| cosmo.inverse_codist(*d))
        .collect()
}

/// Calculate the Rvir from a given mass for a range of redshift values.
/// @param max_solar_mass The maximum viral mass in solar masses for the viral radii.
/// @param redshift_array an array of multiple redshift values.
/// @param omega_m Mass density (often 0.3 in LCDM).
/// @param omega_k Effective mass density of relativistic particles (often 0. in LCDM).
/// @param omega_l Effective mass density of dark energy (often 0.7 in LCDM).
/// @param h0 H0 = 100 * h.
/// @export
#[extendr]
fn calculate_max_rvirs(
    max_solar_mass: f64,
    redshift_array: Vec<f64>,
    omega_m: f64,
    omega_k: f64,
    omega_l: f64,
    h0: f64,
) -> Vec<f64> {
    let cosmo = Cosmology {
        omega_m,
        omega_k,
        omega_l,
        h0,
    };
    redshift_array
        .par_iter()
        .map(|z| cosmo.mvir_to_rvir(max_solar_mass, *z))
        .collect()
}

/// Calculate the Sigma from a given mass for a range of redshift values.
/// @param max_solar_mass The maximum viral mass in solar masses for the viral radii.
/// @param redshift_array an array of multiple redshift values.
/// @param omega_m Mass density (often 0.3 in LCDM).
/// @param omega_k Effective mass density of relativistic particles (often 0. in LCDM).
/// @param omega_l Effective mass density of dark energy (often 0.7 in LCDM).
/// @param h0 H0 = 100 * h.
/// @export
#[extendr]
fn calculate_max_sigmas(
    max_solar_mass: f64,
    redshift_array: Vec<f64>,
    omega_m: f64,
    omega_k: f64,
    omega_l: f64,
    h0: f64,
) -> Vec<f64> {
    let cosmo = Cosmology {
        omega_m,
        omega_k,
        omega_l,
        h0,
    };
    redshift_array
        .par_iter()
        .map(|z| cosmo.mvir_to_sigma(max_solar_mass, *z))
        .collect()
}

/// Distance modulus.
/// @param redshift_array an array of multiple redshift values.
/// @param omega_m Mass density (often 0.3 in LCDM).
/// @param omega_k Effective mass density of relativistic particles (often 0. in LCDM).
/// @param omega_l Effective mass density of dark energy (often 0.7 in LCDM).
/// @param h0 H0 = 100 * h.
/// @returns The distance modulus for the given array of redshifts.
/// @export
#[extendr]
fn distance_modulus(
    redshift_array: Vec<f64>,
    omega_m: f64,
    omega_k: f64,
    omega_l: f64,
    h0: f64,
) -> Vec<f64> {
    let cosmo = Cosmology {
        omega_m,
        omega_k,
        omega_l,
        h0,
    };
    redshift_array
        .par_iter()
        .map(|&z| cosmo.distance_modulus(z))
        .collect()
}

/// differential_comoving_volume
/// @param redshift_array an array of multiple redshift values.
/// @param omega_m Mass density (often 0.3 in LCDM).
/// @param omega_k Effective mass density of relativistic particles (often 0. in LCDM).
/// @param omega_l Effective mass density of dark energy (often 0.7 in LCDM).
/// @param h0 H0 = 100 * h.
/// @returns The distance modulus for the given array of redshifts.
/// @export
fn diff_covol(
    redshift_array: Vec<f64>,
    omega_m: f64,
    omega_k: f64,
    omega_l: f64,
    h0: f64,
) -> Vec<f64> {
    let cosmo = Cosmology {
        omega_m,
        omega_k,
        omega_l,
        h0,
    };
    redshift_array
        .par_iter()
        .map(|&z| cosmo.differential_comoving_distance(z))
        .collect()
}

/// finding the links between all galaxies in a brute force way.
/// @description
/// `fof_links_aaron` will determine all connections between galaxies in a survey and return the pairs.
/// @param ra Array of right ascension values.
/// @param dec Array of declination values.
/// @param comoving_distances Array of comoving distances in Mpc.
/// @param linking_lengths An array of individual scaled linking lengths for each galaxy (ignoring r0 and b0).
/// @param b0 The plane-of-sky constant to be scaled.
/// @param r0 The line-of-sight constant to be scaled.
/// @return A dataframe-like object of tuples which represent the link between galaxies (i, j) if they exist.
/// @export
#[extendr]
fn fof_links_aaron(
    ra_array: Vec<f64>,
    dec_array: Vec<f64>,
    comoving_distances: Vec<f64>,
    linking_lengths_pos: Vec<f64>,
    linking_lengths_los: Vec<f64>,
) -> List {
    let links = ffl1(
        ra_array,
        dec_array,
        comoving_distances,
        linking_lengths_pos,
        linking_lengths_los,
    );
    let i_vec: Vec<usize> = links.iter().map(|(x, _)| *x + 1).collect(); // + 1 for R idx
    let j_vec: Vec<usize> = links.iter().map(|(_, y)| *y + 1).collect();
    list![i = i_vec, j = j_vec]
}

/// finding the links between all galaxies in a brute force way.
/// @description
/// `fof_links_fast` will determine all connections between galaxies in a survey and return the pairs.
/// @param ra Array of right ascension values.
/// @param dec Array of declination values.
/// @param comoving_distances Array of comoving distances in Mpc.
/// @param linking_lengths An array of individual scaled linking lengths for each galaxy (ignoring r0 and b0).
/// @param b0 The plane-of-sky constant to be scaled.
/// @param r0 The line-of-sight constant to be scaled.
/// @return A dataframe-like object of tuples which represent the link between galaxies (i, j) if they exist.
/// @export
#[extendr]
fn fof_links_fast(
    ra_array: Vec<f64>,
    dec_array: Vec<f64>,
    comoving_distances: Vec<f64>,
    linking_lengths_pos: Vec<f64>,
    linking_lengths_los: Vec<f64>,
) -> List {
    let links = find_links(
        ra_array,
        dec_array,
        comoving_distances,
        linking_lengths_pos,
        linking_lengths_los,
    );
    let i_vec: Vec<usize> = links.iter().map(|(x, _)| *x + 1).collect(); // + 1 for R idx
    let j_vec: Vec<usize> = links.iter().map(|(_, y)| *y + 1).collect();
    list![i = i_vec, j = j_vec]
}

/// Creates a group catalog from the given arrays.
/// @description
/// This is the R wrapper for the rust functionality that builds the group catalog. There shouldn't
/// be a need to use this in R as a more R-friendly function should be available.
/// @param ra Array of right asscension values.
/// @param dec Array of declination values.
/// @param redshift Array of redshift values.
/// @param magnitudes Array of apparent magnitude values.
/// @param velocity_errors Array of velocity errors.
/// @param group_id Array of the group ids, where -1 represents galaxies not in a group.
/// @return A named list with group properties.
/// @export
#[extendr]
fn create_group_catalog(
    ra: Vec<f64>,
    dec: Vec<f64>,
    redshift: Vec<f64>,
    absolute_magnitudes: Vec<f64>,
    velocity_errors: Vec<f64>,
    group_ids: Vec<i32>,
    omega_m: f64,
    omega_k: f64,
    omega_l: f64,
    h0: f64,
) -> List {
    let catalog = GroupedGalaxyCatalog {
        ra,
        dec,
        redshift,
        absolute_magnitudes,
        velocity_errors,
        group_ids,
    };
    let cosmo = &Cosmology {
        omega_m,
        omega_k,
        omega_l,
        h0,
    };
    let group_catalog = catalog.calculate_group_properties(cosmo);
    list![
        group_id = group_catalog.ids,
        iter_ra = group_catalog.iter_ras,
        iter_dec = group_catalog.iter_decs,
        iter_redshift = group_catalog.iter_redshifts,
        iter_idx = group_catalog
            .iter_idxs
            .iter()
            .map(|id| id + 1)
            .collect::<Vec<usize>>(), // +1 idx for R
        median_redshift = group_catalog.median_redshifts,
        co_dist = group_catalog.distances,
        r50 = group_catalog.r50s,
        r100 = group_catalog.r100s,
        rsigma = group_catalog.rsigmas,
        multiplicity = group_catalog.multiplicity,
        velocity_dispersion_gap = group_catalog.velocity_dispersion_gap,
        velocity_dispersion_gap_err = group_catalog.velocity_dispersion_gap_err,
        mass_proxy = group_catalog.raw_masses,
        bcg_idxs = group_catalog
            .bcg_idxs
            .iter()
            .map(|id| id + 1)
            .collect::<Vec<usize>>(),
        bcg_ras = group_catalog.bcg_ras,
        bcg_decs = group_catalog.bcg_decs,
        bcg_redshifts = group_catalog.bcg_redshifts,
        center_of_light_ras = group_catalog.col_ras,
        center_of_light_decs = group_catalog.col_decs,
        total_absolute_mag = group_catalog.total_absolute_mags,
        flux_proxies = group_catalog.total_flux_proxies
    ]
}

/// Creates a pair catalog from the given arrays.
/// @description
/// This is the R wrapper for the rust functionality that builds the pair catalog. There shouldn't
/// be a need to use this in R as a more R-friendly function should be available.
/// @param ra Array of right asscension values.
/// @param dec Array of declination values.
/// @param redshift Array of redshift values.
/// @param magnitudes Array of apparent magnitude values.
/// @param group_id Array of the group ids, where -1 represents galaxies not in a group.
/// @return A named list with pair properties.
/// @export
#[extendr]
fn create_pair_catalog(
    ra: Vec<f64>,
    dec: Vec<f64>,
    redshift: Vec<f64>,
    absolute_magnitudes: Vec<f64>,
    group_ids: Vec<i32>,
) -> List {
    let catalog = GroupedGalaxyCatalog {
        ra,
        dec,
        redshift,
        absolute_magnitudes,
        velocity_errors: vec![50.; 1], // dummy variable
        group_ids,
    };

    let pair_catalog = catalog.calculate_pair_properties();
    list![
        pair_id = pair_catalog.ids,
        idx_1 = pair_catalog
            .idx_1
            .iter()
            .map(|id| id + 1)
            .collect::<Vec<i32>>(), // +1 for R
        idx_2 = pair_catalog
            .idx_2
            .iter()
            .map(|id| id + 1)
            .collect::<Vec<i32>>(),
        projected_separation = pair_catalog.projected_separation,
        velocity_separation = pair_catalog.velocity_separation,
        ra_bar = pair_catalog.ra_bar,
        dec_bar = pair_catalog.dec_bar,
        redshift_bar = pair_catalog.redshift_bar,
        total_absolute_mag = pair_catalog.total_absolute_mags,
    ]
}

/// Determines the s score as in Robotham+2011
/// @export
#[extendr]
fn calculate_s_score(measured_groups: &[i32], mock_groups: &[i32], min_group_size: usize) -> f64 {
    s_score(measured_groups, mock_groups, min_group_size)
}

/// The Harmoic Mean
/// @param values The array of objects of which to find the harmonic mean
/// @returns The harmonic mean.
/// @export
#[extendr]
fn calculate_harmonic_mean(values: Vec<f64>) -> f64 {
    harmonic_mean(values)
}

/// Completeness function to calculate the completeness at given positions.
/// @param ra_observed The Right Ascension of the galaxies that were observed, in degrees.
/// @param dec_observed The Declination of the galaxies in degrees that were observed, in degrees.
/// @param ra_target The Right Ascension of the galaxies that were targeted for observation in degrees.
/// @param dec_target The Declination of the galaxies that were targeted for observation in degress.
/// @param angular_radius The search area for each galaxy in degress.
/// @return An array of elements between 0-1.
/// @export
#[extendr]
fn calc_completeness_rust(
    ra_observed: Vec<f64>,
    dec_observed: Vec<f64>,
    ra_target: Vec<f64>,
    dec_target: Vec<f64>,
    angular_radius: Vec<f64>,
) -> Vec<f64> {
    let observed_catalog = PositionCatalog {
        ra_deg: ra_observed,
        dec_deg: dec_observed,
    };
    let target_catalog = PositionCatalog {
        ra_deg: ra_target,
        dec_deg: dec_target,
    };
    calculate_completeness(observed_catalog, target_catalog, angular_radius)
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod Nessie;
    fn comoving_distances_at_z;
    fn z_at_comoving_distances;
    fn distance_modulus;
    fn h_at_z;
    fn calculate_max_rvirs;
    fn calculate_max_sigmas;
    fn fof_links_aaron;
    fn fof_links_fast;
    fn create_group_catalog;
    fn calculate_s_score;
    fn calculate_harmonic_mean;
    fn create_pair_catalog;
    fn calc_completeness_rust;
    fn gen_randoms;
    fn diff_covol;
}

#[cfg(test)]
mod tests {
    use std::iter::zip;

    use super::*;

    #[test]
    fn test_cosdists_basic_lcdm() {
        let z = vec![0.3, 0.5, 0.1, 2.];
        let omega_m = 0.3;
        let omega_k = 0.0;
        let omega_l = 0.7;
        let h0 = 70.;
        let res = comoving_distances_at_z(z, omega_m, omega_k, omega_l, h0);
        let ans = vec![1194.39686972, 1888.62539593, 418.45448763, 5179.86207441];
        for (r, a) in zip(res, ans) {
            assert!((r - a).abs() < 1e-3)
        }
    }

    #[test]
    fn test_inverse_cosdist() {
        let redshifts = vec![0.3, 0.5, 0.1, 2.];
        let omega_m = 0.3;
        let omega_k = 0.0;
        let omega_l = 0.7;
        let h0 = 70.;
        let distances = comoving_distances_at_z(redshifts, omega_m, omega_k, omega_l, h0);
        let returned_z = z_at_comoving_distances(distances, omega_m, omega_k, omega_l, h0);
        let redshifts = vec![0.3, 0.5, 0.1, 2.];
        for (ret_z, z) in zip(returned_z, redshifts) {
            assert!((ret_z - z).abs() < 1e-5)
        }
    }
}
