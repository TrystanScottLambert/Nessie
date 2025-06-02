use extendr_api::prelude::*;
use kiddo::KdTree;
use kiddo::SquaredEuclidean;
use rayon::prelude::*;
use integrate::adaptive_quadrature;
use libm::{tan, atan2};
use std::collections::HashSet;
use std::f64;
use std::iter::zip;
use std::sync::Mutex;
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

// Combined function that finds indices and removes target in one go
fn find_indices_in_range(sorted_array: &[f64], argsort: &[usize], low_lim: f64, up_lim: f64, exclude: usize) -> Vec<usize> {
    
    let start_idx = sorted_array.partition_point(|&x| x < low_lim);
    if start_idx >= sorted_array.len() {
        return Vec::new();
    }
    
    let end_idx = sorted_array.partition_point(|&x| x <= up_lim);
    if start_idx >= end_idx {
        return Vec::new();
    }

    let mut result = Vec::with_capacity(end_idx - start_idx - 1);
    for &idx in &argsort[start_idx..end_idx] {
        if idx > exclude {
            result.push(idx);
        }
    }
    result
}


fn argsort<T: PartialOrd>(data: &[T]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..data.len()).collect();
    idx.sort_by(|&i, &j| data[i].partial_cmp(&data[j]).unwrap());
    idx
}

fn fof_links_rust_old(ra_array: Vec<f64>, dec_array: Vec<f64>, comoving_distances: Vec<f64>, linking_lengths: Vec<f64>, b0: f64, r0: f64) -> Vec<(usize, usize)> {
    
    let n = ra_array.len();
    // we can get away with working out max ll because we precalculate this when working out the
    // linking lengths. Mvir -> Rvir. This is why we don't need to pass a max val and can just calculate it.
    let max_ll = linking_lengths.iter().cloned().fold(f64::NAN, f64::max);
    let mut sorted_distances = comoving_distances.clone();
    sorted_distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let dist_argsort = argsort(&comoving_distances);

    // Option 1: Using Mutex for shared state (simpler but potentially slower)
    let links = Mutex::new(Vec::new());
    
    (0..n).into_par_iter().for_each(|i| {
        let ra_i = ra_array[i];
        let dec_i = dec_array[i];
        let dist_i = comoving_distances[i];
        let ll_i = linking_lengths[i];
        let max_delta_los = r0 * b0 * (max_ll + ll_i) * 0.5;
        let low_lim = dist_i - max_delta_los;
        let up_lim = dist_i + max_delta_los;
        let possible_idx = find_indices_in_range(&sorted_distances, &dist_argsort, low_lim, up_lim, i);

        let mut local_links = Vec::new();
        for j in possible_idx {
            let ra_j = ra_array[j];
            let dec_j = dec_array[j];
            let dist_j = comoving_distances[j];
            let ll_j = linking_lengths[j];
            let avg_link_length = (ll_i + ll_j)/2.;

            if (dist_i - dist_j).abs() <= b0 * r0 * avg_link_length {
                let angular_separation = ang_sep_radians(&ra_i, &dec_i, &ra_j, &dec_j);
                if tan(angular_separation).abs() * (dist_i + dist_j)/2. < b0 * avg_link_length {
                    local_links.push((i, j));
                }
            }
        }
        
        // Only lock once per thread to add all local results
        if !local_links.is_empty() {
            links.lock().unwrap().extend(local_links);
        }
    });
    
    links.into_inner().unwrap()
}



#[derive(Debug, Clone, Copy)]
pub struct SkyCoord {
    pub ra: f64,  // Right Ascension in degrees
    pub dec: f64, // Declination in degrees
}

#[derive(Debug, Clone, Copy)]
pub struct CartesianCoord {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl SkyCoord {
    /// Create a new SkyCoord from RA/DEC in degrees
    pub fn new(ra_deg: f64, dec_deg: f64) -> Self {
        Self {
            ra: ra_deg,
            dec: dec_deg,
        }
    }

    /// Convert to unit sphere Cartesian coordinates (x, y, z)
    pub fn to_cartesian(&self) -> CartesianCoord {
        let ra_rad = self.ra.to_radians();
        let dec_rad = self.dec.to_radians();
        let cos_dec = dec_rad.cos();
        
        CartesianCoord {
            x: cos_dec * ra_rad.cos(),
            y: cos_dec * ra_rad.sin(),
            z: dec_rad.sin(),
        }
    }
}


/// A catalog with a pre-built KD-tree for efficient searching
pub struct SkyCatalog {
    coords: Vec<SkyCoord>,
    tree: KdTree<f64, 3>,
}

impl SkyCatalog {
    /// Create a new catalog from sky coordinates
    pub fn new(ra_array_deg: Vec<f64>, dec_array_deg: Vec<f64>) -> Self {
        let mut coords = Vec::new();
        for (ra, dec) in zip(ra_array_deg, dec_array_deg) {
            coords.push(SkyCoord::new(ra, dec));
        }

        let mut tree = KdTree::new();
        
        for (i, coord) in coords.iter().enumerate() {
            let cart = coord.to_cartesian();
            tree.add(&[cart.x, cart.y, cart.z], i as u64);
        }
        
        Self { coords, tree }
    }

    /// Search for all coordinates within a given angular distance of a specific coordinate
    /// 
    /// # Arguments
    /// * `index` - Index of the coordinate to search around
    /// * `radius_deg` - Search radius in degrees
    /// 
    /// # Returns
    /// * Vector of indices of coordinates within the search radius (excluding the query coordinate itself)
    pub fn search_around(&self, index: usize, radius_deg: f64) -> Vec<usize> {
        let query_coord = &self.coords[index];
        let cart = query_coord.to_cartesian();
        let query_point = [cart.x, cart.y, cart.z];
        
        // Convert angular radius to chord distance
        let radius_rad = radius_deg.to_radians();
        let chord_radius = (2.0 * (radius_rad / 2.0).sin()).powi(2);
        
        // Find all points within the chord distance
        let neighbors = self.tree.within::<SquaredEuclidean>(&query_point, chord_radius);
        // Filter out the query point itself and return just the indices
        neighbors
            .into_iter()
            .map(|neighbor| neighbor.item)
            .filter(|&idx | idx as usize != index)
            .map(|idx| idx as usize)
            .collect()
    }
}


fn fof_links_rust(ra_array: Vec<f64>, dec_array: Vec<f64>, comoving_distances: Vec<f64>, plane_of_sky_linking_lengths: Vec<f64>, line_of_sight_linking_lengths: Vec<f64>) -> Vec<(usize, usize)> {
    
    let n = ra_array.len();
    // we can get away with working out max ll because we precalculate this when working out the
    // linking lengths. Mvir -> Rvir. This is why we don't need to pass a max val and can just calculate it.
    let max_pos_ll = plane_of_sky_linking_lengths.iter().cloned().fold(f64::NAN, f64::max);
    let max_los_ll = line_of_sight_linking_lengths.iter().cloned().fold(f64::NAN, f64::max);
    let mut sorted_distances = comoving_distances.clone();
    sorted_distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let dist_argsort = argsort(&comoving_distances);
    let sky_catalog = SkyCatalog::new(ra_array.clone(), dec_array.clone());

    // Option 1: Using Mutex for shared state (simpler but potentially slower)
    let links = Mutex::new(Vec::new());
    
    (0..n).into_par_iter().for_each(|i| {
        let ra_i = ra_array[i];
        let dec_i = dec_array[i];
        let dist_i = comoving_distances[i];

        let maximum_angle_radians = ((2. * max_pos_ll)/(dist_i + sorted_distances[0])).atan();

        let low_lim = dist_i - max_los_ll;
        let up_lim = dist_i + max_los_ll;
        let possible_los_idx = find_indices_in_range(&sorted_distances, &dist_argsort, low_lim, up_lim, i);
        
        let possible_pos_idx: HashSet<usize> = sky_catalog.search_around(i, maximum_angle_radians.to_degrees())
            .into_iter()
            .filter(|idx| idx > &i )
            .collect();

        let possible_idx: Vec<usize> = possible_los_idx
            .into_iter()
            .filter(|x| possible_pos_idx
            .contains(x))
            .collect();

        let mut local_links = Vec::new();
        for j in possible_idx {
            let ra_j = ra_array[j];
            let dec_j = dec_array[j];
            let dist_j = comoving_distances[j];

            if (dist_i - dist_j).abs() <= 0.5 * (line_of_sight_linking_lengths[i] + line_of_sight_linking_lengths[j]) {
                let angular_separation = ang_sep_radians(&ra_i, &dec_i, &ra_j, &dec_j);
                if tan(angular_separation).abs() * (dist_i + dist_j)/2. < 0.5 * (plane_of_sky_linking_lengths[i] + plane_of_sky_linking_lengths[j]) {
                    local_links.push((i, j));
                }
            }
        }
        
        // Only lock once per thread to add all local results
        if !local_links.is_empty() {
            links.lock().unwrap().extend(local_links);
        }
    });
    
    links.into_inner().unwrap()
}

fn ang_sep_radians(lon_1: &f64, lat_1: &f64, lon_2: &f64, lat_2: &f64) -> f64 {
    // we should do a global check that all ra and dec and z are well behaved.
    // function is assuming nothing dumb is happening. 
    let lon1 = lon_1.to_radians();
    let lat1 = lat_1.to_radians();
    let lon2 = lon_2.to_radians();
    let lat2 = lat_2.to_radians();

    let sdlon = (lon2 - lon1).sin();
    let cdlon = (lon2 - lon1).cos();
    let slat1 = (lat1).sin();
    let slat2 = (lat2).sin();
    let clat1 = (lat1).cos();
    let clat2 = (lat2).cos();

    let num1 = clat2 * sdlon;
    let num2 = clat1 * slat2 - slat1 * clat2 * cdlon;
    let denominator = slat1 * slat2 + clat1 * clat2 * cdlon;
    let hypot = (num1.powi(2) + num2.powi(2)).sqrt();
    atan2(hypot, denominator)
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
fn fof_links(ra_array: Vec<f64>, dec_array: Vec<f64>, comoving_distances: Vec<f64>, linking_lengths_pos: Vec<f64>, linking_lengths_los: Vec<f64>) -> List {
    let links = fof_links_rust(ra_array, dec_array, comoving_distances, linking_lengths_pos, linking_lengths_los);
    let i_vec: Vec<usize> = links.iter().map(|(x, _)| *x + 1).collect(); // + 1 for R idx
    let j_vec: Vec<usize> = links.iter().map(|(_, y)|  *y + 1).collect();
    list![i = i_vec, j = j_vec]
}


fn fof_links_brutal(ra_array: Vec<f64>, dec_array: Vec<f64>, comoving_distances: Vec<f64>, linking_lengths_pos: Vec<f64>, linking_lengths_los: Vec<f64>) -> Vec<(usize, usize)> {

    let mut links: Vec<(usize, usize)> = Vec::new();
    let n = ra_array.len();
    for i in 0..(n - 1) {
        for j in (i + 1)..(n) {
            let avg_link_length = (linking_lengths_los[i] + linking_lengths_los[j])/2.;
            if (comoving_distances[i] - comoving_distances[j]).abs() <= avg_link_length {
                let angular_separation = ang_sep_radians(&ra_array[i], &dec_array[i], &ra_array[j], &dec_array[j]);
                if tan(angular_separation).abs() * (comoving_distances[i] + comoving_distances[j])/2. < (linking_lengths_pos[i] + linking_lengths_pos[j]) * 0.5 {
                    links.push((i, j));
                }
            }
        }
    }
    links
}

/// finding the links between all galaxies in a brute force way.
/// @description
/// `fof_link_brutal` will determine all connections between galaxies in a survey and return the pairs.
/// @param ra Array of right ascension values.
/// @param dec Array of declination values.
/// @param comoving_distances Array of comoving distances in Mpc. 
/// @param linking_lengths An array of individual scaled linking lengths for each galaxy (ignoring r0 and b0).
/// @param b0 The plane-of-sky constant to be scaled. 
/// @param r0 The line-of-sight constant to be scaled.
/// @return A dataframe-like object of tuples which represent the link between galaxies (i, j) if they exist.
/// @export
#[extendr]
fn fof_links_brute_force(ra_array: Vec<f64>, dec_array: Vec<f64>, comoving_distances: Vec<f64>, linking_lengths_pos: Vec<f64>, linking_lengths_los: Vec<f64>) -> List {
    let links = fof_links_brutal(ra_array, dec_array, comoving_distances, linking_lengths_pos, linking_lengths_los);
    let i_vec: Vec<usize> = links.iter().map(|(x, _)| *x + 1).collect(); // + 1 for R idx
    let j_vec: Vec<usize> = links.iter().map(|(_, y)|  *y + 1).collect();
    list![i = i_vec, j = j_vec]
}


pub fn ffl1(
    ra_array: Vec<f64>,
    dec_array: Vec<f64>,
    comoving_distances: Vec<f64>,
    linking_lengths_pos: Vec<f64>,
    linking_lengths_los: Vec<f64>,
) -> Vec<(usize, usize)> {
    let n = ra_array.len();

    // Convert (RA, Dec, dist) to 3D Cartesian coordinates
    let coords: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let ra_rad = ra_array[i].to_radians();
            let dec_rad = dec_array[i].to_radians();
            let x = dec_rad.cos() * ra_rad.cos();
            let y = dec_rad.cos() * ra_rad.sin();
            let z = dec_rad.sin();
            [x, y, z]
        })
        .collect();
    
    let mut ind = Vec::new();

    for i in 0..(n - 1) {
        for j in (i+1)..n {
            let ztemp = (linking_lengths_los[i] + linking_lengths_los[j]) * 0.5;
            let zrad = (comoving_distances[i] - comoving_distances[j]).abs();

            if zrad <= ztemp {
                let bgal2 = ((linking_lengths_pos[i] + linking_lengths_pos[j]) * 0.5).powi(2);

                let radproj = (0..3)
                    .map(|k| (coords[i][k] - coords[j][k]).powi(2))
                    .sum::<f64>();

                if radproj <= bgal2 {
                    ind.push((i, j));
                }
            }
        }
    }
    ind
}

/// finding the links between all galaxies in a brute force way.
/// @description
/// `fof_link_brutal` will determine all connections between galaxies in a survey and return the pairs.
/// @param ra Array of right ascension values.
/// @param dec Array of declination values.
/// @param comoving_distances Array of comoving distances in Mpc. 
/// @param linking_lengths An array of individual scaled linking lengths for each galaxy (ignoring r0 and b0).
/// @param b0 The plane-of-sky constant to be scaled. 
/// @param r0 The line-of-sight constant to be scaled.
/// @return A dataframe-like object of tuples which represent the link between galaxies (i, j) if they exist.
/// @export
#[extendr]
fn fof_links_aaron(ra_array: Vec<f64>, dec_array: Vec<f64>, comoving_distances: Vec<f64>, linking_lengths_pos: Vec<f64>, linking_lengths_los: Vec<f64>) -> List {
    let links = ffl1(ra_array, dec_array, comoving_distances, linking_lengths_pos, linking_lengths_los);
    let i_vec: Vec<usize> = links.iter().map(|(x, _)| *x + 1).collect(); // + 1 for R idx
    let j_vec: Vec<usize> = links.iter().map(|(_, y)| *y + 1).collect();
    list![i = i_vec, j = j_vec]
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod Nessie;
    fn comoving_distances_at_z;
    fn z_at_comoving_distances;
    fn fof_links;
    fn fof_links_brute_force;
    fn fof_links_aaron;
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
            seps.push(ang_sep_radians(ra1, dec1, ra2, dec2).to_degrees());
        }

        for (a_sep, sep) in zip(seps_astropy, seps) {
            assert!((a_sep - sep).abs() < 1e-5)
        }
    }

}