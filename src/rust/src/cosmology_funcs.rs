
use integrate::adaptive_quadrature;
use libm::log10;
use libm::sinh;
use roots::SimpleConvergency;
use roots::find_root_brent;
use rayon::prelude::*;


use crate::constants::SPEED_OF_LIGHT;

pub fn e_func(z: f64, omega_m: f64, omega_k: f64, omega_l: f64) -> f64 {
    (omega_m * (1.0 + z).powi(3) + omega_k * (1.0 + z).powi(2) + omega_l).sqrt()
}


pub fn hubble_distance(hubble_constant: f64) -> f64 {
    SPEED_OF_LIGHT/hubble_constant
}


pub fn comoving_distance(redshift: f64, omega_m: f64, omega_k: f64, omega_l: f64, h0: f64) -> f64 {
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


pub fn inverse_codist(distance: f64, omega_m: f64, omega_k: f64, omega_l: f64, h0: f64) -> f64 {
    let f = |z: f64| {comoving_distance(z, omega_m, omega_k, omega_l, h0) - distance};
    let mut convergency = SimpleConvergency {eps:1e-8f64, max_iter: 30};
    match find_root_brent(1e-9, 1200., &f, &mut convergency) {
        Ok(t) => t,
        Err(_error) => 0.0
    }
}

#[allow(dead_code)]
pub fn comoving_transverse_distance(redshift: f64, omega_m: f64, omega_k: f64, omega_l: f64, h0: f64) -> f64 {
    let co_dist = comoving_distance(redshift, omega_m, omega_k, omega_l, h0);
    let h_dist = hubble_distance(h0);

    match omega_k {
        val if val > 0. => {h_dist * (1./omega_k.sqrt()) * sinh(omega_k.sqrt() * (co_dist/h_dist))},
        val if val < 0. => {h_dist * (1./omega_k.abs().sqrt()) * (omega_k.abs().sqrt() * (co_dist/h_dist)).sin()},
        _ => co_dist

    }
}

pub fn distance_modulus(redshift: f64, omega_m: f64, omega_k: f64, omega_l: f64, h0: f64) -> f64 {
    let co_trans_dist = comoving_transverse_distance(redshift, omega_m, omega_k, omega_l, h0);
    5. * log10(co_trans_dist * (1. + redshift)) + 25.
}

pub fn dist_mods(redshifts: Vec<f64>, omega_m: f64, omega_k: f64, omega_l: f64, h0: f64) -> Vec<f64> {
    redshifts
        .par_iter()
        .map(|z| distance_modulus(*z, omega_m, omega_l, omega_k, h0))
        .collect()
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_e_func_basic_lcdm() {
        let z = 0.3;
        let omega_m = 0.3;
        let omega_k = 0.0;
        let omega_l = 0.7;
        let e = e_func(z, omega_m, omega_k, omega_l);
        assert!((e - 1.16580444329227).abs() < 1e-5);
    }

}