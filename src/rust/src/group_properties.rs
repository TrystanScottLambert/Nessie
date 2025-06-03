use std::{f64::consts::PI, iter::zip};
use stats::{mean, median};

use crate::spherical_trig_funcs::{convert_equitorial_to_cartesian,convert_cartesian_to_equitorial, convert_equitorial_to_cartesian_scaled, euclidean_distance_3d};
use crate::constants::SPEED_OF_LIGHT;
use crate::cosmology_funcs::{distance_modulus};
use crate::helper_funcs::quantile_interpolated;


// sigma error would be in km/s. 
pub fn velocity_dispersion_gapper(group_redshifts: Vec<f64>, group_redshifts_err: Vec<f64>) -> (f64, f64) {
    let sigma_err_squared = mean(group_redshifts_err.iter().copied());
    let median_redshift = median(group_redshifts.iter().copied()).unwrap();
    let n  = group_redshifts.len();
    let nf64 = n as f64;

    let mut velocities: Vec<f64> = group_redshifts
        .iter()
        .map(|z| (z*SPEED_OF_LIGHT)/(1. + median_redshift))
        .collect();

    velocities.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let gaps: Vec<f64> = velocities.windows(2).map(|w| w[1] - w[0]).collect();
    let i: Vec<usize> = (1..n).collect();
    let weights: Vec<usize> = i.iter().map(|val| val * (n - val)).collect();

    let sum: f64 = zip(weights, gaps).map(|(a, b)| a as f64 * b).sum();
    let sigma_gap: f64 = ((PI.sqrt())/(nf64 * (nf64 - 1.))) * sum;
    let raw_dispersion_squared = (n as f64 * sigma_gap.powi(2)) / (nf64 - 1.);
    let dispersion = (raw_dispersion_squared - sigma_err_squared).sqrt();
    (dispersion, sigma_err_squared.sqrt())
}


pub fn calculate_iterative_center(ra_group: Vec<f64>, dec_group: Vec<f64>, redshift_group: Vec<f64>, magnitudes_group: Vec<f64>,  omega_m: f64, omega_k: f64, omega_l: f64, h0: f64) -> (f64, f64) {
    
    let absolute_magnitudes: Vec<f64> = zip(magnitudes_group, redshift_group)
        .map(|(mag, z)| mag - distance_modulus(z, omega_m, omega_k, omega_l, h0))
        .collect();
    let coords_cartesian: Vec<[f64; 3]> = zip(ra_group, dec_group)
        .map(|(ra, dec)| convert_equitorial_to_cartesian(&ra, &dec))
        .collect();
    let flux: Vec<f64> = absolute_magnitudes
        .iter()
        .map(|mag| 10_f64.powf(-0.4 * mag))
        .collect();

    let mut temp_flux = flux.clone();
    let mut temp_coords = coords_cartesian.clone();

    while temp_flux.len() > 2 {
        let flux_sum: f64 = temp_flux.iter().cloned().sum();

        let center: [f64; 3] = (0..3).map(|i| {
            temp_coords.iter().zip(temp_flux.iter())
                .map(|(coord, &f)| coord[i] * f).sum::<f64>() / flux_sum
        }).collect::<Vec<f64>>().try_into().unwrap();

        let distances: Vec<f64> = temp_coords.iter()
            .map(|coord| {
            ((coord[0] - center[0]).powi(2) + 
            (coord[1] - center[1]).powi(2) +
            (coord[2] - center[2]).powi(2)).sqrt()
            }).collect();

        if let Some((max_idx, _)) = distances.iter().enumerate().max_by(|a, b| a.1.partial_cmp(b.1).unwrap()) {
            temp_coords.remove(max_idx);
            temp_flux.remove(max_idx);
        } else {
            break;
        }
    }

    let max_flux_idx = temp_flux
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    let final_cartesian = temp_coords[max_flux_idx];

    // Convert back to spherical RA/Dec (in degrees)
    let center = convert_cartesian_to_equitorial(&final_cartesian[0], &final_cartesian[1], &final_cartesian[2]);
    let wrapped_ra = if center[0] < 0.0 { center[0] + 360.0 } else { center[0] };

    (wrapped_ra, center[1])
}


pub fn calculate_radius(
    group_ra: Vec<f64>,
    group_dec: Vec<f64>,
    group_dists: Vec<f64>,
    group_center_ra: f64,
    group_center_dec: f64,
    group_center_dist: f64,
) -> [f64; 3] {
    let center = convert_equitorial_to_cartesian_scaled(group_center_ra, group_center_dec, group_center_dist);

    let mut distances: Vec<f64> = group_ra
        .iter()
        .zip(&group_dec)
        .zip(&group_dists)
        .map(|((&ra, &dec), &d)| {
            let pos = convert_equitorial_to_cartesian_scaled(ra, dec, d);
            euclidean_distance_3d(&pos, &center)
        })
        .collect();

    distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
    println!("Distances {:?}", distances);

    let q50 = quantile_interpolated(&distances, 0.5);
    let q68 = quantile_interpolated(&distances, 0.68);
    let q100 = *distances.last().unwrap();

    [q50, q68, q100]
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gapper_velocity_basic() {
        let redshifts = vec![0.3, 0.2, 0.32, 0.5];
        let err = vec![50., 40., 30., 10.];
        let result = velocity_dispersion_gapper(redshifts, err);
        assert_eq!(result.0, 35908.750028805);
        assert_eq!(result.1, 5.70087712549569);
    }


    #[test]
    fn test_bright_central_galaxy_dominates_center() {
        let ra_group = vec![180.0, 179.0, 181.0, 180.0];
        let dec_group = vec![0.0, 1.0, -1.0, 0.0];
        let redshift_group = vec![0.1; 4];
        let magnitudes_group = vec![12.0, 18.0, 18.0, 18.0]; // central galaxy is much brighter

        let omega_m = 0.3;
        let omega_k = 0.0;
        let omega_l = 0.7;
        let h0 = 70.0;

        let (ra, dec) = calculate_iterative_center(
            ra_group,
            dec_group,
            redshift_group,
            magnitudes_group,
            omega_m,
            omega_k,
            omega_l,
            h0
        );

        // RA should be close to 180.0 and Dec to 0.0 (within ~0.01 deg)
        assert!((ra - 180.0).abs() < 1e-6, "RA deviated too far: {}", ra);
        assert!(dec.abs() < 1e-6, "Dec deviated too far: {}", dec);
    }


    #[test]
    fn comparing_radii_to_r() {
        let group_ra = 23.;
        let group_dec = -23.;
        let group_dist = 20.;

        let group_ras = vec![23., 23.2, 22.9, 24.0];
        let group_decs = vec![-23., -23.2, -23.2, -23.];
        let group_dists = vec![20., 20., 20., 20.];

        let answer =  [0.08584888, 0.10391348, 0.32131273];
        let result = calculate_radius(group_ras, group_decs, group_dists, group_ra, group_dec, group_dist);
        for (r, a) in zip(result, answer) {
            assert!((r - a).abs() < 1e-6)
        }
    }
}