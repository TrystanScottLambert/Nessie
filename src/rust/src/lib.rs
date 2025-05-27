use extendr_api::prelude::*;
use rayon::prelude::*;
use integrate::adaptive_quadrature;
use libm::{tan, atan2};
use std::f64::{self, consts::PI};
use roots::SimpleConvergency;
use roots::find_root_brent;


const SPEED_OF_LIGHT: f64 = 299_792.458; // km/s

/// Evaluates the E(z) function often used in determining other cosmological functions.
/// @param z single redshift value.
/// @param omega_m Mass density (often 0.3 in LCDM)
/// @param omega_k Effective mass density of relativistic particles (often 0. in LCDM)
/// @param omega_l Effective mass density of dark energy (often 0.7 in LCDM)
/// @return E(z)
/// @export
#[extendr]
fn e_func(z: f64, omega_m: f64, omega_k: f64, omega_l: f64) -> f64 {
    (omega_m * (1.0 + z).powi(3) + omega_k * (1.0 + z).powi(2) + omega_l).sqrt()
}

/// Hubble time multiplied by the speed of light.
/// @param hubble_constant H0 = 100 * h
/// @return c/hubble_constant
/// @export
#[extendr]
fn hubble_distance(hubble_constant: f64) -> f64 {
    SPEED_OF_LIGHT/hubble_constant
}

/// Calculates the comoving distance for a single redshift value.
/// @param z single redshift value.
/// @param omega_m Mass density (often 0.3 in LCDM)
/// @param omega_k Effective mass density of relativistic particles (often 0. in LCDM)
/// @param omega_l Effective mass density of dark energy (often 0.7 in LCDM)
/// @param hubble_constant H0 = 100 * h
/// @return a single comoving distance in Mpc.
/// @export
#[extendr]
fn comoving_distance(redshift: f64, omega_m: f64, omega_k: f64, omega_l: f64, h0: f64) -> f64 {
    if redshift < 1e-7 {
        return 0.;
    }
    let tolerance = 1e-5;
    let min_h = 1e-7;
    let f = |z:f64| 1./e_func(z, omega_m, omega_k, omega_l);
    let cosmo_recession_velocity = adaptive_quadrature::adaptive_simpson_method(f, 0.0, redshift, min_h, tolerance)
        .expect("Value too close to zero. Must be within 10e-8");
    hubble_distance(h0) * cosmo_recession_velocity
}

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

/// redshift at some given comoving distance in Mpc
/// @export
#[extendr]
fn inverse_codist(distance: f64, omega_m: f64, omega_k: f64, omega_l: f64, h0: f64) -> f64 {
    let f = |z: f64| {comoving_distance(z, omega_m, omega_k, omega_l, h0) - distance};
    let mut convergency = SimpleConvergency {eps:1e-8f64, max_iter: 30};
    match find_root_brent(1e-9, 1200., &f, &mut convergency) {
        Ok(t) => t,
        Err(_error) => 0.0
    }
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


fn ang_sep(lon_1: &f64, lat_1: &f64, lon_2: &f64, lat_2: &f64) -> f64 {
    assert!(*lon_1 >= 0., "lon_1 must be in [0, 360)");
    assert!(*lon_2 >= 0., "lon_2 must be in [0, 360)");
    assert!((-90.0..=90.0).contains(lat_1), "lat_1 must be in [-90, 90]");
    assert!((-90.0..=90.0).contains(lat_2), "lat_2 must be in [-90, 90]");

    let sdlon = (lon_2.to_radians() - lon_1.to_radians()).sin();
    let cdlon = (lon_2.to_radians() - lon_1.to_radians()).cos();
    let slat1 = (lat_1.to_radians()).sin();
    let slat2 = (lat_2.to_radians()).sin();
    let clat1 = (lat_1.to_radians()).cos();
    let clat2 = (lat_2.to_radians()).cos();

    let num1 = clat2 * sdlon;
    let num2 = clat1 * slat2 - slat1 * clat2 * cdlon;
    let denominator = slat1 * slat2 + clat1 * clat2 * cdlon;
    let hypot = (num1.powi(2) + num2.powi(2)).sqrt();
    atan2(hypot, denominator).to_degrees()
}

/// finding the links between all galaxies
/// @description
/// `fof_links` will determine all connections between galaxies in a survey and return the pairs.
/// @param ra Array of right ascension values.
/// @param dec Array of declination values.
/// @param comoving_distances Array of comoving distances in Mpc. 
/// @param linking_lengths An array of individual scaled linking lengths for each galaxy (ignoring r0 and b0).
/// @param b0 The plane-of-sky constant to be scaled. 
/// @param r0 The line-of-sight constant to be scaled.
/// @return A dataframe-like object of tuples which represent the link between galaxies (i, j) if they exist.
/// @export
#[extendr]
fn fof_links(ra_array: &[f64], dec_array: &[f64], comoving_distances: &[f64], linking_lengths: &[f64], b0: f64, r0: f64) -> List {

    let mut link_i: Vec<usize> = Vec::new();
    let mut link_j: Vec<usize> = Vec::new();

    for i in 0..(ra_array.len() - 1) {
        for j in (i + 1)..(ra_array.len()) {
            let avg_link_length = (linking_lengths[i] + linking_lengths[j])/2.;
            if (comoving_distances[i] - comoving_distances[j]).abs() <= b0 * r0 * avg_link_length {
                let angular_separation = ang_sep(&ra_array[i], &dec_array[i], &ra_array[j], &dec_array[j]);
                if angular_separation.tan() * (comoving_distances[i] + comoving_distances[j])/2. < b0 * avg_link_length {
                    link_i.push(i+1); // R indexing
                    link_j.push(j+1);
                }
            }
        }
    }
    list![i = link_i, j = link_j]
}



// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod Nessie;
    fn comoving_distances_at_z;
    fn z_at_comoving_distances;
    fn fof_links;
}


#[cfg(test)]
mod tests {

    use super::*;
    use std::iter::zip;

    #[test]
    fn test_e_func_basic_lcdm() {
        let z = 0.3;
        let omega_m = 0.3;
        let omega_k = 0.0;
        let omega_l = 0.7;
        let e = e_func(z, omega_m, omega_k, omega_l);
        assert!((e - 1.16580444329227).abs() < 1e-5);
    }

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
            print!("{ret_z} {z}");
            assert!((ret_z - z).abs() < 1e-5)
        }
    }

    #[test]
    fn test_angsep() {
        let ras_1 = [0., 90., 270., 360., 0.3];
        let ras_2 = [12., 23., 44., 180., 0.2];
        let decs_1 = [-90., 12., -23., 0., 1.];
        let decs_2 = [-88., 4., 4., -2., 80.];

        let seps_astropy = vec![2.,  66.68631254, 131.69266712, 178.,79.00001543];
        let mut seps: Vec<f64> = Vec::new();
        
        for i in 0..ras_1.len() {
            let ra1 = ras_1.get(i).unwrap();
            let ra2 = ras_2.get(i).unwrap();
            let dec1 = decs_1.get(i).unwrap();
            let dec2 = decs_2.get(i).unwrap();
            seps.push(ang_sep(ra1, dec1, ra2, dec2));
        }

        for (a_sep, sep) in zip(seps_astropy, seps) {
            assert!((a_sep - sep).abs() < 1e-5)
        }
    }

    #[test]
    fn test_fof_linking_onsky_condition() {
        let ras = [121.1, 121.2, 121.3, 181.1, 181.2, 181.1];
        let decs = [-23.1, -23.1, -23.1, 68., 68., 68.3];
        let distances = [20., 20., 20., 20., 20., 20.];
        let linking_lengths = [2., 2., 2., 2., 2., 2.];
        // with the above settings b0 = tan(0.2deg)*10 = 0.0349 
        // should connect [(0:2), (3:4), (5)]
        let result = fof_links(&ras, &decs, &distances, &linking_lengths, 0.0349, 123.);
        
    }
}