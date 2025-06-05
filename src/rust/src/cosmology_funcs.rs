
use integrate::adaptive_quadrature;
use libm::log10;
use libm::sinh;
use roots::SimpleConvergency;
use roots::find_root_brent;

use crate::constants::SPEED_OF_LIGHT;


pub struct Cosmology {
    pub omega_m: f64,
    pub omega_k: f64,
    pub omega_l: f64,
    pub h0: f64
}

impl Cosmology {
    pub fn e_func(&self, z: f64) -> f64 {
        (self.omega_m * (1.0 + z).powi(3) + self.omega_k * (1.0 + z).powi(2) + self.omega_l).sqrt()
    }

    pub fn hubble_distance(&self) -> f64 {
        SPEED_OF_LIGHT/self.h0
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
        

}



#[cfg(test)]
mod tests {
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

        println!("{}",cosmo.distance_modulus(z));
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

}