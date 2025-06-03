
use libm::{atan2, asin};

// assuming a unit sphere
pub fn convert_equitorial_to_cartesian(ra_deg: &f64, dec_deg: &f64) -> [f64; 3] {
    let ra_radians = ra_deg.to_radians();
    let dec_radians = dec_deg.to_radians();
    let x = dec_radians.cos() * ra_radians.cos();
    let y = dec_radians.cos() * ra_radians.sin();
    let z = dec_radians.sin();
    [x, y, z]
}

// assuming a unit sphere
pub fn convert_cartesian_to_equitorial(x: &f64, y: &f64, z: &f64) -> [f64; 2] {
    let radius = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
    let ra_radian = atan2(*y, *x);
    let dec_radian = asin(z/radius);
    [ra_radian.to_degrees(), dec_radian.to_degrees()]
}

