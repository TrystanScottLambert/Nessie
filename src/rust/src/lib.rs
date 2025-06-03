use extendr_api::prelude::*;
use rayon::prelude::*;

use std::f64;

pub mod group_properties;
pub mod spherical_trig_funcs;
pub mod cosmology_funcs;
pub mod constants;
use crate::cosmology_funcs::{comoving_distance, inverse_codist};
use crate::spherical_trig_funcs::{convert_equitorial_to_cartesian};


/// Calculates multiple comoving distances for multiple redshifts.
/// @param redshift_array an array of multiple redshift values.
/// @param omega_m Mass density (often 0.3 in LCDM)
/// @param omega_k Effective mass density of relativistic particles (often 0. in LCDM)
/// @param omega_l Effective mass density of dark energy (often 0.7 in LCDM)
/// @param hubble_constant H0 = 100 * h
/// @return multiple comoving distance in Mpc.
/// @export
#[extendr]
fn comoving_distances_at_z(redshift_array: Vec<f64>, omega_m:f64, omega_k:f64, omega_l:f64, h0: f64) -> Vec<f64> {
    redshift_array
        .par_iter()
        .map(|z| comoving_distance(*z, omega_m, omega_k, omega_l, h0))
        .collect()
}


/// redshift at some given comoving distances in Mpc
/// @export
#[extendr]
fn z_at_comoving_distances(distances: Vec<f64>, omega_m: f64, omega_k: f64, omega_l: f64, h0: f64) -> Vec<f64> {
    distances
        .par_iter()
        .map(|d| inverse_codist(*d, omega_m, omega_k, omega_l, h0))
        .collect()
}


pub fn ffl1(
    ra_array: Vec<f64>,
    dec_array: Vec<f64>,
    comoving_distances: Vec<f64>,
    linking_lengths_pos: Vec<f64>,
    linking_lengths_los: Vec<f64>,
) -> Vec<(usize, usize)> {
    let n = ra_array.len();

    // Convert (RA, Dec, dist) to 3D Cartesian coordinates
    let coords: Vec<[f64; 3]> = (0..n)
        .map(|i| convert_equitorial_to_cartesian(&ra_array[i], &dec_array[i]))
        .collect();
    
    let mut ind = Vec::new();

    for i in 0..(n - 1) {
        for j in (i+1)..n {
            let ztemp = (linking_lengths_los[i] + linking_lengths_los[j]) * 0.5;
            let zrad = (comoving_distances[i] - comoving_distances[j]).abs();

            if zrad <= ztemp {
                let bgal2 = ((linking_lengths_pos[i] + linking_lengths_pos[j]) * 0.5).powi(2);

                let radproj = (0..3)
                    .map(|k| (coords[i][k] - coords[j][k]).powi(2))
                    .sum::<f64>();

                if radproj <= bgal2 {
                    ind.push((i, j));
                }
            }
        }
    }
    ind
}

/// finding the links between all galaxies in a brute force way.
/// @description
/// `fof_link_brutal` will determine all connections between galaxies in a survey and return the pairs.
/// @param ra Array of right ascension values.
/// @param dec Array of declination values.
/// @param comoving_distances Array of comoving distances in Mpc. 
/// @param linking_lengths An array of individual scaled linking lengths for each galaxy (ignoring r0 and b0).
/// @param b0 The plane-of-sky constant to be scaled. 
/// @param r0 The line-of-sight constant to be scaled.
/// @return A dataframe-like object of tuples which represent the link between galaxies (i, j) if they exist.
/// @export
#[extendr]
fn fof_links_aaron(ra_array: Vec<f64>, dec_array: Vec<f64>, comoving_distances: Vec<f64>, linking_lengths_pos: Vec<f64>, linking_lengths_los: Vec<f64>) -> List {
    let links = ffl1(ra_array, dec_array, comoving_distances, linking_lengths_pos, linking_lengths_los);
    let i_vec: Vec<usize> = links.iter().map(|(x, _)| *x + 1).collect(); // + 1 for R idx
    let j_vec: Vec<usize> = links.iter().map(|(_, y)| *y + 1).collect();
    list![i = i_vec, j = j_vec]
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod Nessie;
    fn comoving_distances_at_z;
    fn z_at_comoving_distances;
    fn fof_links_aaron;
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