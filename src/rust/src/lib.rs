use extendr_api::prelude::*;
use ndarray::Array2;
use rayon::prelude::*;
use integrate::adaptive_quadrature;
use libm::{sin, sinh, log10, asinh};
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
    if redshift == 0. {
        return 0.;
    }
    let tolerance = 10.0e-6;
    let min_h = 10.0e-8;
    let f = |z:f64| 1./e_func(z, omega_m, omega_k, omega_l);
    let cosmo_recession_velocity = adaptive_quadrature::adaptive_simpson_method(f, 0.0, redshift, min_h, tolerance)
        .expect("Value to close to zero. Must be within 10e-8");
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



#[derive(Debug)]
struct FoFResult {
    pairs: Vec<(usize, usize)>,
    radial_distances: Vec<f64>,
    redshift_differences: Vec<f64>,
    truncated: bool,
}



fn fof_linking(
    positions: &Array2<f64>,
    bgal: &[f64],
    zcut: &[f64],
    nmax: usize,
) -> FoFResult {
    let n = positions.nrows();

    let mut pairs = Vec::with_capacity(nmax.min(n * (n - 1) / 2));
    let mut radial_distances = Vec::with_capacity(nmax.min(n * (n - 1) / 2));
    let mut redshift_differences = Vec::with_capacity(nmax.min(n * (n - 1) / 2));

    for i in 0..n {
        for j in (i + 1)..n {
            let z_thresh = 0.5 * (zcut[i] + zcut[j]);
            let z_diff = (positions[[i, 3]] - positions[[j, 3]]).abs();

            if z_diff > z_thresh {
                continue;
            }

            let b2 = 0.5 * (bgal[i] + bgal[j]);
            let b2 = b2 * b2;

            let r2 = (0..3)
                .map(|k| (positions[[i, k]] - positions[[j, k]]).powi(2))
                .sum::<f64>();

            if r2 <= b2 {
                if pairs.len() >= nmax {
                    return FoFResult {
                        pairs,
                        radial_distances,
                        redshift_differences,
                        truncated: true,
                    };
                }

                pairs.push((i + 1, j + 1)); // 1-based indexing
                radial_distances.push(r2.sqrt());
                redshift_differences.push(z_diff);
            }
        }
    }

    FoFResult {
        pairs,
        radial_distances,
        redshift_differences,
        truncated: false,
    }
}


// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod Nessie;
    fn comoving_distances_at_z;
    fn z_at_comoving_distances;
}