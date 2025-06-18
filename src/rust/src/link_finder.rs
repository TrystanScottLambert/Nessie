use rayon::prelude::*;

use crate::spherical_trig_funcs::convert_equitorial_to_cartesian;



// Combined function that finds indices and removes target in one go
fn find_indices_in_range(sorted_array: &[f64], argsort_array: &[usize], low_lim: f64, up_lim: f64, exclude_idx: usize) -> Vec<usize> {
    
    let start_idx = sorted_array.partition_point(|&x| x < low_lim);
    if start_idx >= sorted_array.len() {
        return Vec::new();
    }
    
    let end_idx = sorted_array.partition_point(|&x| x <= up_lim);
    if start_idx >= end_idx {
        return Vec::new();
    }

    let mut result = Vec::with_capacity(end_idx - start_idx - 1);
    for &idx in &argsort_array[start_idx..end_idx] {
        if idx > exclude_idx {
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
        .map(|i| convert_equitorial_to_cartesian(&ra_array[i], &dec_array[i]))
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

pub fn fast_ffl1(
    ra_array: Vec<f64>,
    dec_array: Vec<f64>,
    comoving_distances: Vec<f64>,
    linking_lengths_pos: Vec<f64>,
    linking_lengths_los: Vec<f64>,
) -> Vec<(usize, usize)> {
    let n = ra_array.len();

    // Convert (RA, Dec, dist) to 3D Cartesian coordinates
    let coords: Vec<[f64; 3]> = (0..n)
        .map(|i| convert_equitorial_to_cartesian(&ra_array[i], &dec_array[i]))
        .collect();
    
    let max_los_ll = linking_lengths_los.iter().cloned().fold(f64::NAN, f64::max);
    let mut sorted_distances = comoving_distances.clone();
    sorted_distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let dist_argsort = argsort(&comoving_distances);

    let mut ind = Vec::new();

    for i in 0..(n - 1) {
        let dist_i = comoving_distances[i];
        let lower_lim = dist_i - max_los_ll;
        let upper_lim = dist_i + max_los_ll;
        let possible_los_idx = find_indices_in_range(&sorted_distances, &dist_argsort, lower_lim, upper_lim, i);

        for j in possible_los_idx {
            let average_los_ll = (linking_lengths_los[i] + linking_lengths_los[j]) * 0.5;
            let zrad = (comoving_distances[i] - comoving_distances[j]).abs();

            if zrad <= average_los_ll {
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





pub fn fast_ffl1_with_spatial_grid(
    ra_array: Vec<f64>,
    dec_array: Vec<f64>,
    comoving_distances: Vec<f64>,
    linking_lengths_pos: Vec<f64>,
    linking_lengths_los: Vec<f64>,
) -> Vec<(usize, usize)> {
    let n = ra_array.len();
    
    // Convert coordinates
    let coords: Vec<[f64; 3]> = (0..n)
        .map(|i| convert_equitorial_to_cartesian(&ra_array[i], &dec_array[i]))
        .collect();

    // Create spatial grid
    let max_pos_ll = linking_lengths_pos.iter().cloned().fold(f64::NAN, f64::max);
    let grid_size = max_pos_ll * 2.0; // Adjust based on your data
    
    // Build spatial hash map
    let mut spatial_grid: std::collections::HashMap<(i32, i32, i32), Vec<usize>> = 
        std::collections::HashMap::new();
    
    for (i, coord) in coords.iter().enumerate() {
        let grid_x = (coord[0] / grid_size).floor() as i32;
        let grid_y = (coord[1] / grid_size).floor() as i32;
        let grid_z = (coord[2] / grid_size).floor() as i32;
        
        spatial_grid.entry((grid_x, grid_y, grid_z))
            .or_default()
            .push(i);
    }

    // Sort by comoving distance
    let mut sorted_indices: Vec<usize> = (0..n).collect();
    sorted_indices.sort_unstable_by(|&i, &j| {
        comoving_distances[i].partial_cmp(&comoving_distances[j]).unwrap()
    });

    let max_los_ll = linking_lengths_los.iter().cloned().fold(f64::NAN, f64::max);

    // Process in parallel
    let results: Vec<(usize, usize)> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            let mut local_results = Vec::new();
            
            // Get grid cell for particle i
            let grid_x = (coords[i][0] / grid_size).floor() as i32;
            let grid_y = (coords[i][1] / grid_size).floor() as i32;
            let grid_z = (coords[i][2] / grid_size).floor() as i32;
            
            // Check neighboring grid cells
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let key = (grid_x + dx, grid_y + dy, grid_z + dz);
                        if let Some(neighbors) = spatial_grid.get(&key) {
                            for &j in neighbors {
                                if j == i {
                                    continue;
                                }
                                
                                // Quick comoving distance pre-filter using max_los_ll
                                let zrad = (comoving_distances[i] - comoving_distances[j]).abs();
                                if zrad > max_los_ll {
                                    continue;
                                }
                                
                                // Precise line-of-sight check
                                let los_ll = 0.5 * (linking_lengths_los[i] + linking_lengths_los[j]);
                                if zrad > los_ll {
                                    continue;
                                }
                                
                                // Spatial distance check
                                let radproj2 = (coords[i][0] - coords[j][0]).powi(2)
                                    + (coords[i][1] - coords[j][1]).powi(2)
                                    + (coords[i][2] - coords[j][2]).powi(2);
                                
                                let bgal2 = 0.5 * (linking_lengths_pos[i] + linking_lengths_pos[j]);
                                let bgal2 = bgal2 * bgal2;
                                
                                if radproj2 <= bgal2 {
                                    local_results.push((i, j));
                                }
                            }
                        }
                    }
                }
            }
            
            local_results
        })
        .collect();

    results
}



pub fn fast_ffl1_parallel(
    ra_array: Vec<f64>,
    dec_array: Vec<f64>,
    comoving_distances: Vec<f64>,
    linking_lengths_pos: Vec<f64>,
    linking_lengths_los: Vec<f64>,
) -> Vec<(usize, usize)> {
    let n = ra_array.len();

    let coords: Vec<[f64; 3]> = (0..n)
        .map(|i| convert_equitorial_to_cartesian(&ra_array[i], &dec_array[i]))
        .collect();

    let max_los_ll = linking_lengths_los
        .iter()
        .cloned()
        .fold(f64::NAN, f64::max);

    let mut sorted_distances = comoving_distances.clone();
    sorted_distances.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let dist_argsort = argsort(&comoving_distances);

    // Parallel outer loop
    let results: Vec<(usize, usize)> = (0..(n - 1))
        .into_par_iter()
        .map(|i| {
            let dist_i = comoving_distances[i];
            let lower_lim = dist_i - max_los_ll;
            let upper_lim = dist_i + max_los_ll;
            let possible_los_idx =
                find_indices_in_range(&sorted_distances, &dist_argsort, lower_lim, upper_lim, i);

            let mut local_pairs = Vec::new();

            for j in possible_los_idx {
                let average_los_ll = (linking_lengths_los[i] + linking_lengths_los[j]) * 0.5;
                let zrad = (comoving_distances[i] - comoving_distances[j]).abs();

                if zrad <= average_los_ll {
                    let bgal2 = ((linking_lengths_pos[i] + linking_lengths_pos[j]) * 0.5).powi(2);

                    let radproj = (0..3)
                        .map(|k| (coords[i][k] - coords[j][k]).powi(2))
                        .sum::<f64>();

                    if radproj <= bgal2 {
                        local_pairs.push((i, j));
                    }
                }
            }

            local_pairs
        })
        .flatten()
        .collect();

    results
}




#[cfg(test)]
mod tests {
    use std::iter::zip;

    use super::*;

    #[test]
    fn faster_links() {
        // These functions should be better benchmarked in the R code reading in the 
        // old gama data.

        // Doesn't really matter what we put here. They have to match FFL1 exactly.
        // 3 2 1 group set up. 
        let ra_array = vec![20., 20., 20., -20., -20., 0.];
        let dec_array = vec![-50., -50., -50., 50., 50., 75.];
        let comoving_distances = vec![100., 100., 100., 20., 20., 200.];
        let linking_lengths_pos = vec![0.06; 6];
        let linking_lengths_los = vec![18.; 6];

        let classic_links = ffl1(ra_array.clone(), dec_array.clone(), comoving_distances.clone(), linking_lengths_pos.clone(), linking_lengths_los.clone());
        let new_links = fast_ffl1_parallel(ra_array.clone(), dec_array.clone(), comoving_distances.clone(), linking_lengths_pos.clone(), linking_lengths_los.clone());
        for (new, classic) in zip(new_links, classic_links) {
            assert_eq!(new.0, classic.0);
            assert_eq!(new.1, classic.1);
        }
    }


    #[test]
    fn argsort_simple() {
        let data = [5., 3., 4., 1., 3.];
        let answer: [usize; 5] = [3, 1, 4, 2, 0];
        let result = argsort(&data);
        for (a, r) in zip(answer, result) {
            assert_eq!(a, r)
        }
    }

    #[test]
    fn test_find_indices_in_range() {
        let _data = [5., 3., 4., 1., 3., 0., 10.];
        let ans = [4];
        let data_sorted = [0., 1., 3., 3., 4., 5., 10.];
        let argsort = [5, 3, 1, 4, 2, 0, 6];
        let lower_lim = 3.;
        let upper_lim = 6.;
        
        let idx = find_indices_in_range(&data_sorted, &argsort, lower_lim, upper_lim, 2);
        assert_eq!(idx, &ans);
    }



}