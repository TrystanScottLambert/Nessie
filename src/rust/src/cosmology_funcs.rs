
use std::f64::consts::PI;

use integrate::adaptive_quadrature;
use libm::log10;
use libm::sinh;
use roots::SimpleConvergency;
use roots::find_root_brent;

use crate::constants::SPEED_OF_LIGHT;

pub const G: f64 = 6.67384e-11; // m^3 kg^-1 s^-2
pub const KM_TO_METERS: f64 = 1000.;
pub const MPC_TO_METERS: f64 = 3.08567758e22;
pub const MSOL_TO_KG: f64 = 1.9891e30;
pub const PC_TO_METERS: f64 = 3.08568e16;
  

pub struct Cosmology {
    pub omega_m: f64,
    pub omega_k: f64,
    pub omega_l: f64,
    pub h0: f64
}

#[derive(PartialEq)]
pub enum DistanceType {
    Angular,
    Comoving
}

impl Cosmology {
    pub fn e_func(&self, z: f64) -> f64 {
        (self.omega_m * (1.0 + z).powi(3) + self.omega_k * (1.0 + z).powi(2) + self.omega_l).sqrt()
    }

    pub fn hubble_distance(&self) -> f64 {
        SPEED_OF_LIGHT/self.h0
    }

    pub fn h_at_z(&self, z: f64) -> f64 {
        self.h0 * self.e_func(z)
    }
    
    pub fn rho_critical(&self, z: f64, distance_type: DistanceType) -> f64 {
        let hub_square = self.h_at_z(z).powi(2);
        let mut rho_crit = ((3. * hub_square) / (8. * PI * G)) * (KM_TO_METERS.powi(2) / MPC_TO_METERS.powi(2));
        
        rho_crit /= MSOL_TO_KG;
        rho_crit *= MPC_TO_METERS.powi(3);


        match distance_type {
            DistanceType::Angular => rho_crit,
            DistanceType::Comoving => rho_crit / (1. + z).powi(3)
        }
        
    }

    pub fn comoving_distance(&self, z: f64) -> f64 {
        if z < 1e-7 {
            return 0.;
        }
        let tolerance = 1e-5;
        let min_h = 1e-7;
        let f = |z:f64| 1./self.e_func(z);
        let cosmo_recession_velocity = adaptive_quadrature::adaptive_simpson_method(f, 0.0, z, min_h, tolerance)
            .expect("Value too close to zero. Must be within 10e-8");
        self.hubble_distance() * cosmo_recession_velocity
    }


    pub fn inverse_codist(&self, distance: f64) -> f64 {
        let f = |z: f64| {self.comoving_distance(z) - distance};
        let mut convergency = SimpleConvergency {eps:1e-8f64, max_iter: 30};
        match find_root_brent(1e-9, 1200., &f, &mut convergency) {
            Ok(t) => t,
            Err(_error) => 0.0
        }
    }

    pub fn comoving_transverse_distance(&self, z: f64) -> f64 {
        let co_dist = self.comoving_distance(z);
        let h_dist = self.hubble_distance();

        match self.omega_k {
            val if val > 0. => {h_dist * (1./self.omega_k.sqrt()) * sinh(self.omega_k.sqrt() * (co_dist/h_dist))},
            val if val < 0. => {h_dist * (1./self.omega_k.abs().sqrt()) * (self.omega_k.abs().sqrt() * (co_dist/h_dist)).sin()},
            _ => co_dist

        }
    }

    pub fn distance_modulus(&self, z: f64) -> f64 {
        let co_trans_dist = self.comoving_transverse_distance(z);
        5. * log10(co_trans_dist * (1. + z)) + 25.
    }

    pub fn mvir_to_rvir(&self, solar_mass: f64, z: f64) -> f64 {
        let rho_crit = self.rho_critical(z, DistanceType::Comoving); // 1e9/1e6 for Lunit and other units. 
        (solar_mass * (3./4.) * (1./(PI * 200. * rho_crit))).powf(1./3.) // 200 for delta_vir = 200
    }

    pub fn mvir_to_sigma(&self, solar_mass: f64, z: f64) -> f64 {
        let rho_crit = self.rho_critical(z, DistanceType::Angular);
        let g = G * (MSOL_TO_KG) / (PC_TO_METERS);
        let unit_constant = g/1e12;
        (solar_mass * (32.*PI*(unit_constant.powi(3))* 200. * rho_crit/3.).sqrt()).powf(1./3.) // 200 for delta_vir = 200
    }

}


#[cfg(test)]
mod tests {
    use std::iter::zip;

    use super::*;

    #[test]
    fn test_e_func_basic_lcdm() {
        let z = 0.3;
        let cosmo = Cosmology {
            omega_m: 0.3,
            omega_k: 0.,
            omega_l: 0.7,
            h0: 70.,
        };

        let e = cosmo.e_func(z);
        assert!((e - 1.16580444329227).abs() < 1e-5);
    }

    #[test]
    fn testing_flat_cosmo_versus_celestial() {
        let z = 0.3;
        let cosmo = Cosmology {
            omega_m: 0.3,
            omega_k: 0.,
            omega_l: 0.7,
            h0: 70.,
        };

        assert!((cosmo.comoving_distance(z) - 1194.397).abs() < 1e-3);
        assert!((cosmo.comoving_transverse_distance(z) - 1194.397).abs() < 1e-3);
        assert!((cosmo.distance_modulus(z) - 40.95546).abs() < 1e-3);
    }

    #[test]
    fn testing_inverse_codist_function() {
        let z = 0.3;
        let cosmo = Cosmology {
            omega_m: 0.3,
            omega_k: 0.,
            omega_l: 0.7,
            h0: 70.,
        };

        let co_dist = cosmo.comoving_distance(z);
        let inverse = cosmo.inverse_codist(co_dist);
        assert!((inverse - z).abs() < 1e-7);
    }

    #[test]
    fn testing_h_at_z() {
        // comparing to celestials h grow function with a flat cosmology.
        let z = 0.3;
        let cosmo = Cosmology {
            omega_m: 0.3,
            omega_k: 0.,
            omega_l: 0.7,
            h0: 100.,
        };

        let result = cosmo.h_at_z(z);
        let answer = 116.5804;
        assert!((result - answer).abs() < 1e-4);

    }

    #[test]
    fn testing_rho_critical() {
        // Comparing to celestials cosgrowRhoCrit function for a flat cosmology with Dist = 'Ang'
        let z = 0.3;
        let cosmo = Cosmology {
            omega_m: 0.3,
            omega_k: 0.,
            omega_l: 0.7,
            h0: 100.,
        };

        // angular
        let answer = 377095148067.;
        let result = cosmo.rho_critical(z, DistanceType::Angular).trunc();
        assert_eq!(answer, result);

        // co moving
        let answer = 171640941314.;
        let result = cosmo.rho_critical(z, DistanceType::Comoving).trunc();
        assert_eq!(answer, result);

    }

    #[test]
    fn testing_mvir_to_rvir() {
        // Comparing to celestial coshaloMvirtoRvir
        let z = 0.3;
        let cosmo = Cosmology {
            omega_m: 0.3,
            omega_k: 0.,
            omega_l: 0.7,
            h0: 100.,
        };
        let solar_mass = 1e12;
        let result = cosmo.mvir_to_rvir(solar_mass, z);
        let answer = 0.190877;
        assert!((result - answer).abs() < 1e-5);

    }

    #[test]
    fn testing_mvir_to_sigma() {
        // Comparing to celestial coshaloMvirtoRvir
        let z = 0.3;
        let cosmo = Cosmology {
            omega_m: 0.3,
            omega_k: 0.,
            omega_l: 0.7,
            h0: 100.,
        };
        let solar_mass = 1e12;
        let result = cosmo.mvir_to_sigma(solar_mass, z);
        let answer = 242.07541;
        assert!((result - answer).abs() < 1e-5);
    }

    #[test]
    fn testing_distance_modulus() {
        // Comparing this calculation to celestials distance modulus.
        let cosmo = Cosmology {omega_m: 0.3, omega_k: 0., omega_l: 0.7, h0:100.};
        let redshifts = [0.1, 0.2, 0.3, 1., 2., 4.];
        let answers = [37.54069, 39.18177, 40.18095, 43.32573, 45.18269, 46.99805];
        let results: Vec<f64> = redshifts.iter().map(|&z| cosmo.distance_modulus(z)).collect();
        for (r, a) in zip(results, answers) {
            assert!((r - a).abs() < 1e-5)
        }
    }


}