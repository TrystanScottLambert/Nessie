use std::{f64::consts::PI, iter::zip};
use stats::{mean, median};

use crate::spherical_trig_funcs::{convert_equitorial_to_cartesian,convert_cartesian_to_equitorial, convert_equitorial_to_cartesian_scaled, euclidean_distance_3d};
use crate::constants::SPEED_OF_LIGHT;
use crate::cosmology_funcs::Cosmology;
use crate::helper_funcs::quantile_interpolated;

pub struct Group {
    pub ra_members: Vec<f64>,
    pub dec_members: Vec<f64>,
    pub redshift_members: Vec<f64>,
    pub magnitudes_members: Vec<f64>,
    pub velocity_errors: Vec<f64>
}

impl Group {
    // sigma error would be in km/s. 
    pub fn velocity_dispersion_gapper(&self) -> (f64, f64) {
        let sigma_err_squared = mean(self.velocity_errors.iter().copied());
        let median_redshift = median(self.redshift_members.iter().copied()).unwrap();
        let n  = self.redshift_members.len();
        let nf64 = n as f64;

        let mut velocities: Vec<f64> = self.redshift_members
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



    pub fn calculate_iterative_center(&self, cosmo: Cosmology) -> (f64, f64) {
    
        let absolute_magnitudes: Vec<f64> = zip(self.magnitudes_members.clone(), self.redshift_members.clone())
            .map(|(mag, z)| mag - cosmo.distance_modulus(z))
            .collect();
        let coords_cartesian: Vec<[f64; 3]> = zip(self.ra_members.clone(), self.dec_members.clone())
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



    pub fn calculate_radius(&self, group_center_ra: f64, group_center_dec: f64, group_center_z: f64, cosmo: Cosmology) -> [f64; 3] {
        let group_center_dist= cosmo.comoving_distance(group_center_z);
        let group_dists: Vec<f64> = self.redshift_members.iter().map(|z| cosmo.comoving_distance(*z)).collect();
        let center = convert_equitorial_to_cartesian_scaled(group_center_ra, group_center_dec, group_center_dist);

        let mut distances: Vec<f64> = self.ra_members
            .iter()
            .zip(&self.dec_members)
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

}



#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_gapper_velocity_basic() {

        let group = Group {
            redshift_members: vec![0.3, 0.2, 0.32, 0.5],
            velocity_errors: vec![50., 40., 30., 10.],
            ra_members: vec![0.2, 0.2, 0.2, 0.2],
            dec_members: vec![50., 50., 50., 50.],
            magnitudes_members: vec![18., 18., 18., 18.]
        };


        let result = group.velocity_dispersion_gapper();
        assert_eq!(result.0, 35908.750028805);
        assert_eq!(result.1, 5.70087712549569);
    }


    #[test]
    fn test_bright_central_galaxy_dominates_center() {
        let group = Group{
            ra_members: vec![180.0, 179.0, 181.0, 180.0],
            dec_members: vec![0.0, 1.0, -1.0, 0.0],
            redshift_members: vec![0.1; 4],
            magnitudes_members: vec![12.0, 18.0, 18.0, 18.0], // central very bright.
            velocity_errors: vec![50., 50., 50., 50.]
        };

        let cosmo = Cosmology {
            omega_m: 0.3,
            omega_k: 0.,
            omega_l: 0.7,
            h0: 0.7
        };

        let (ra, dec) = group.calculate_iterative_center(cosmo);
        // RA should be close to 180.0 and Dec to 0.0 (within ~0.01 deg)
        assert!((ra - 180.0).abs() < 1e-6, "RA deviated too far: {}", ra);
        assert!(dec.abs() < 1e-6, "Dec deviated too far: {}", dec);
    }


    #[test]
    fn comparing_radii_to_r() {
        let cosmo = Cosmology {
            omega_m: 0.3,
            omega_k: 0.,
            omega_l: 0.7,
            h0: 0.7
        };

        let group = Group {
            ra_members: vec![23., 23.2, 22.9, 24.0],
            dec_members: vec![-23., -23.2, -23.2, -23.],
            redshift_members: vec![0.2, 0.2, 0.2, 0.2],
            magnitudes_members: vec![18., 18., 18., 18.],
            velocity_errors: vec![50., 50., 50., 50.]
        };

        let group_ra = 23.;
        let group_dec = -23.;
        let group_z = 0.2;

        let answer =  [350.5738627578198, 424.34275593189506, 1312.1178128946078];
        let result = group.calculate_radius(group_ra, group_dec, group_z, cosmo);
        for (r, a) in zip(result, answer) {
            assert!((r - a).abs() < 1e-6)
        }
    }
}