
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

// Converts RA, Dec (degrees) and comoving distance to 3D Cartesian coordinates.
pub fn convert_equitorial_to_cartesian_scaled(ra_deg: f64, dec_deg: f64, distance: f64) -> [f64; 3] {
    let ra_rad = ra_deg.to_radians();
    let dec_rad = dec_deg.to_radians();
    let x = distance * dec_rad.cos() * ra_rad.cos();
    let y = distance * dec_rad.cos() * ra_rad.sin();
    let z = distance * dec_rad.sin();
    [x, y, z]
}

// The three dimensional euclidean distance
pub fn euclidean_distance_3d(point_1: &[f64; 3], point_2: &[f64; 3]) -> f64 {
    ((point_1[0] - point_2[0]).powi(2) +
    (point_1[1] - point_2[1]).powi(2) +
    (point_1[2] - point_2[2]).powi(2)).sqrt()
}