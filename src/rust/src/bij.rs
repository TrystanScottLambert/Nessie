use std::collections::HashMap;
use rayon::prelude::*;

#[derive(Debug)]
pub struct BijResults {
    pub e_num: usize,
    pub e_den: usize,
    pub q_num: f64,
    pub q_den: f64,
}

pub fn bijcheck(group_ids_1: &[i32], group_ids_2: &[i32], min_group_size: usize) -> BijResults {
    assert_eq!(group_ids_1.len(), group_ids_2.len());

    // Frequency tables excluding -1 - these need to be computed sequentially
    let mut count_table_1 = HashMap::<i32, usize>::new();
    let mut count_table_2 = HashMap::<i32, usize>::new();

    for &g in group_ids_1.iter().filter(|&&x| x != -1) {
        *count_table_1.entry(g).or_insert(0) += 1;
    }
    for &g in group_ids_2.iter().filter(|&&x| x != -1) {
        *count_table_2.entry(g).or_insert(0) += 1;
    }

    // Filter groups in tab1 with size >= groupcut
    let valid_groups_1: Vec<i32> = count_table_1
        .iter()
        .filter_map(|(&group, &count)| if count >= min_group_size { Some(group) } else { None })
        .collect();

    // Indices of valid group members - can be parallelized
    let valid_indices_1: Vec<usize> = group_ids_1
        .par_iter()
        .enumerate()
        .filter(|(_, &g)| valid_groups_1.contains(&g))
        .map(|(i, _)| i)
        .collect();

    // Create group_list sequentially to maintain order
    let group_list: Vec<i32> = {
        let mut seen = HashMap::new();
        valid_indices_1
            .iter()
            .filter_map(|&i| {
                let g = group_ids_1[i];
                if seen.insert(g, ()).is_none() {
                    Some(g)
                } else {
                    None
                }
            })
            .collect()
    };

    // Parallel processing of each group
    let results: Vec<(f64, f64, f64)> = group_list
        .par_iter()
        .map(|&group_id| {
            // Find all galaxies in this group
            let group_galaxies: Vec<usize> = group_ids_1
                .iter()
                .enumerate()
                .filter_map(|(idx, &g)| if g == group_id { Some(idx) } else { None })
                .collect();

            let overlap_groups: Vec<i32> = group_galaxies.iter().map(|&idx| group_ids_2[idx]).collect();
            let overlap_valid: Vec<i32> = overlap_groups.iter().cloned().filter(|&g| g != -1).collect();

            let n1_current = *count_table_1.get(&group_id).unwrap_or(&1) as f64;

            let (q1, q2) = if !overlap_valid.is_empty() {
                let mut temptab = HashMap::<i32, usize>::new();
                for g in &overlap_valid {
                    *temptab.entry(*g).or_insert(0) += 1;
                }

                let mut frac_1 = Vec::new();
                let mut frac_2 = Vec::new();

                for (&group2, &count) in &temptab {
                    if let Some(&n2_val) = count_table_2.get(&group2) {
                        frac_1.push(count as f64 / n1_current);
                        frac_2.push(count as f64 / n2_val as f64);
                    }
                }

                let num_isolated = overlap_groups.iter().filter(|&&g| g == -1).count();
                if num_isolated > 0 {
                    let iso_frac1 = 1.0 / n1_current;
                    for _ in 0..num_isolated {
                        frac_1.push(iso_frac1);
                        frac_2.push(1.0);
                    }
                }

                // Fixed: Use R's which.max() behavior - return FIRST index of maximum value
                // Single pass to find first occurrence of maximum
                let mut best_match = 0;
                let mut max_product = f64::NEG_INFINITY;
                
                for (idx, (&f1, &f2)) in frac_1.iter().zip(&frac_2).enumerate() {
                    let product = f1 * f2;
                    if product > max_product {
                        max_product = product;
                        best_match = idx;
                    }
                    // Note: we don't update best_match if product == max_product,
                    // which gives us the "first occurrence" behavior
                }

                (frac_1[best_match], frac_2[best_match])
            } else {
                // All isolated
                (1.0 / n1_current, 1.0)
            };

            (q1, q2, n1_current)
        })
        .collect();

    // Extract results back into separate vectors
    let mut q1_values = Vec::with_capacity(results.len());
    let mut q2_values = Vec::with_capacity(results.len());
    let mut n1_values = Vec::with_capacity(results.len());

    for (q1, q2, n1) in results {
        q1_values.push(q1);
        q2_values.push(q2);
        n1_values.push(n1);
    }

    let e_num = q1_values
        .iter()
        .zip(&q2_values)
        .filter(|(&q1, &q2)| q1 > 0.5 && q2 > 0.5)
        .count();

    let e_den = n1_values.len();
    let q_num: f64 = q1_values
        .iter()
        .zip(&n1_values)
        .map(|(&q1, &n1)| q1 * n1)
        .sum();

    let q_den: f64 = n1_values.iter().sum();
    BijResults {
        e_num,
        e_den,
        q_num,
        q_den,
    }
}

pub fn s_score(measured_groups: &[i32], mock_groups: &[i32], groupcut: usize) -> f64 {
    let mock_vals = bijcheck(mock_groups, measured_groups, groupcut);
    let measured_vals = bijcheck(measured_groups, mock_groups, groupcut);
    let mock_e = mock_vals.e_num as f64 / mock_vals.e_den as f64;
    let fof_e = measured_vals.e_num as f64 / measured_vals.e_den as f64;
    let mock_q = mock_vals.q_num / mock_vals.q_den;
    let fof_q = measured_vals.q_num / measured_vals.q_den;
    println!("mock_e: {}, fof_e: {}, mock_q: {}, fof_q: {}", mock_vals.e_num, measured_vals.e_num, mock_vals.q_num, measured_vals.q_num);
    println!("mock_e_den: {}, fof_e_den: {}, mock_q_den: {}, fof_q_den: {}", mock_vals.e_den, measured_vals.e_den, mock_vals.q_den, measured_vals.q_den);
    mock_e * fof_e * mock_q * fof_q
}