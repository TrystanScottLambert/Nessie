use std::{f64::consts::PI, iter::zip};
use stats::{mean, median};

pub const SPEED_OF_LIGHT: f64 = 299_792.458; // km/s

// sigma error would be in km/s. 
#[allow(dead_code)]
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
}