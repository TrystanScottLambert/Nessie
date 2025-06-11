use std::{collections::{HashMap, HashSet}};
use rayon::prelude::*;


fn count_groups(group_array: &[i32], singleton_id: i32) -> usize {
    let unique_ids: HashSet<i32> = group_array.iter()
        .filter(|&idx|  *idx != singleton_id)
        .copied()
        .collect();
    unique_ids.len()
}

fn count_group_galaxies(group_array: &[i32], singleton_id: i32) -> usize {
    let group_galaxies: Vec<i32> = group_array
        .iter()
        .filter(|&x| *x != singleton_id)
        .cloned()
        .collect();
    group_galaxies.len()
} 


pub struct QuantityResult {
    pub fof_metric: f64,
    pub mock_metric: f64,
    pub total_metric: f64
}

#[derive(Debug, Clone)]
pub struct GroupCalcs {
    pub purity: f64,
    pub length: usize
}

#[derive(Debug)]
pub struct CombinedResults {
    pub bijective_count: usize,
    pub purity_per_group: HashMap<i32, GroupCalcs>,
}


fn calculate_bijective_and_purity(
    fof_groups: &[i32], 
    mock_groups: &[i32], 
    singleton_id: i32
) -> CombinedResults {
    calculate_bijective_and_purity_with_min_size(fof_groups, mock_groups, singleton_id, 2)
}

fn calculate_bijective_and_purity_with_min_size(
    fof_groups: &[i32], 
    mock_groups: &[i32], 
    singleton_id: i32,
    min_group_size: usize,
) -> CombinedResults {

    // Build index maps with minimum group size filtering
    let mut group1_to_indices: HashMap<i32, Vec<usize>> = HashMap::new();
    let mut group2_to_indices: HashMap<i32, Vec<usize>> = HashMap::new();
    let mut temp_group1: HashMap<i32, Vec<usize>> = HashMap::new();
    let mut temp_group2: HashMap<i32, Vec<usize>> = HashMap::new();

    // First pass: collect all groups (including small ones for mock groups)
    for (i, (&g1, &g2)) in fof_groups.iter().zip(mock_groups.iter()).enumerate() {
        if g1 != singleton_id {
            temp_group1.entry(g1).or_default().push(i);
        }
        if g2 != singleton_id {
            temp_group2.entry(g2).or_default().push(i);
        }
    }

    // Second pass: only keep FOF groups with min_group_size+ members
    // But keep ALL mock groups (we need them for overlap calculations)
    for (group_id, indices) in temp_group1 {
        if indices.len() >= min_group_size {
            group1_to_indices.insert(group_id, indices);
        }
    }
    
    // Keep all mock groups, regardless of size
    group2_to_indices = temp_group2;

    
    // Process groups in parallel
    let results: Vec<_> = group1_to_indices
        .par_iter()
        .map(|(&fof_group_id, indices_fof)| {
            let fof_len = indices_fof.len();
            
            // Count all galaxies and track their mock group assignments
            let mut mock_group_counts: HashMap<i32, usize> = HashMap::new();
            let mut num_isolated = 0;
            let mut num_small_groups = 0;
            
            // First pass: count where each galaxy in this FOF group maps to
            for &idx in indices_fof {
                let mock_id = mock_groups[idx];
                if mock_id == singleton_id {
                    num_isolated += 1;
                } else if let Some(mock_indices) = group2_to_indices.get(&mock_id) {
                    if mock_indices.len() >= min_group_size {
                        *mock_group_counts.entry(mock_id).or_insert(0) += 1;
                    } else {
                        // Galaxy belongs to a mock group that's too small
                        num_small_groups += 1;
                    }
                }
            }
            
            let mut best_q1 = 0.0;
            let mut best_q2 = 0.0;
            let mut found_bijective_match = false;
            
            // Calculate fractions for each mock group overlap
            let mut all_fractions: Vec<(f64, f64)> = Vec::new();
            
            for (&mock_id, &overlap_count) in &mock_group_counts {
                if let Some(mock_indices) = group2_to_indices.get(&mock_id) {
                    let mock_len = mock_indices.len();
                    
                    let frac1 = overlap_count as f64 / fof_len as f64;  // Q1: fraction of FOF group
                    let frac2 = overlap_count as f64 / mock_len as f64; // Q2: fraction of mock group
                    
                    all_fractions.push((frac1, frac2));
                }
            }
            
            // Add contributions from isolated galaxies (like the R code padding)
            if num_isolated > 0 {
                let isolated_frac1 = 1.0 / fof_len as f64;  // Each isolated galaxy contributes 1/N1
                let isolated_frac2 = 1.0;                   // Each is 100% of its own "group"
                for _ in 0..num_isolated {
                    all_fractions.push((isolated_frac1, isolated_frac2));
                }
            }
            
            // Add contributions from galaxies in small mock groups (treat as isolated)
            if num_small_groups > 0 {
                let small_group_frac1 = 1.0 / fof_len as f64;
                let small_group_frac2 = 1.0;
                for _ in 0..num_small_groups {
                    all_fractions.push((small_group_frac1, small_group_frac2));
                }
            }
            
            // Find the best match (highest combined score)
            if !all_fractions.is_empty() {
                let mut best_combined = 0.0;
                let mut best_idx = 0;
                
                for (i, &(frac1, frac2)) in all_fractions.iter().enumerate() {
                    let combined = frac1 * frac2;
                    if combined > best_combined {
                        best_combined = combined;
                        best_idx = i;
                    }
                }
                
                best_q1 = all_fractions[best_idx].0;
                best_q2 = all_fractions[best_idx].1;
                
                // Check for bijective match
                if best_q1 > 0.5 && best_q2 > 0.5 {
                    found_bijective_match = true;
                }
            }

            (fof_group_id, GroupCalcs {
                purity: best_q1,  // Use Q1 as the purity score (like G1int contribution)
                length: fof_len,
            }, found_bijective_match, best_q1, best_q2)
        })
        .collect();


    // Collect results
    let mut purity_results: HashMap<i32, GroupCalcs> = HashMap::with_capacity(results.len());
    let mut bijective_matches = 0;
    let mut g1int_num = 0.0;
    let mut g1int_den = 0.0;

    for (fof_group_id, group_calcs, found_bijective, q1, _q2) in results {
        purity_results.insert(fof_group_id, group_calcs.clone());
        
        if found_bijective {
            bijective_matches += 1;
        }
        
        // Calculate G1int components (like the R code)
        g1int_num += q1 * group_calcs.length as f64;
        g1int_den += group_calcs.length as f64;
    }

    CombinedResults {
        bijective_count: bijective_matches,
        purity_per_group: purity_results,
    }
}

pub struct STotalMetrics {
    pub e_metrics: QuantityResult,
    pub q_metrics: QuantityResult,
    pub s_total: f64
}


// Helper function to calculate Q score from purity results
fn calculate_q_score_from_purity(
    purity_results: &HashMap<i32, GroupCalcs>, 
    total_galaxies: usize
) -> f64 {
    let sum: f64 = purity_results
        .values()
        .map(|calcs| (calcs.length as f64) * calcs.purity)
        .sum();
    sum / total_galaxies as f64
}

// Updated function that does everything in one go
pub fn calculate_s_total(fof_groups: &[i32], mock_groups: &[i32], singleton_id: i32) -> STotalMetrics {
    
    // Get both directions in one combined call each
    let fof_to_mock = calculate_bijective_and_purity(fof_groups, mock_groups, singleton_id);
    let mock_to_fof = calculate_bijective_and_purity(mock_groups, fof_groups, singleton_id);
    
    
    // Calculate E metrics
    let n_bij = fof_to_mock.bijective_count as f64;
    let n_fof_groups = count_groups(fof_groups, singleton_id) as f64;
    let n_mock_groups = count_groups(mock_groups, singleton_id) as f64;
    
    let e_fof = n_bij / n_fof_groups;
    let e_mock = n_bij / n_mock_groups;
    let e_metrics = QuantityResult {
        fof_metric: e_fof,
        mock_metric: e_mock,
        total_metric: e_fof * e_mock,
    };
    
    // Calculate Q metrics
    let total_fof_galaxies = count_group_galaxies(fof_groups, singleton_id);
    let total_mock_galaxies = count_group_galaxies(mock_groups, singleton_id);
    
    let q_fof = calculate_q_score_from_purity(&fof_to_mock.purity_per_group, total_fof_galaxies);
    let q_mock = calculate_q_score_from_purity(&mock_to_fof.purity_per_group, total_mock_galaxies);
    let q_metrics = QuantityResult {
        fof_metric: q_fof,
        mock_metric: q_mock,
        total_metric: q_fof * q_mock,
    };
    
    // Calculate S total
    let s_total = e_metrics.total_metric * q_metrics.total_metric;
    
    
    STotalMetrics {
        e_metrics,
        q_metrics,
        s_total,
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_counting_groups() {
        let groups = [1, 1, 1, -1, 2, 2, 3, 3, 3];
        let actual = count_groups(&groups, -1);
        assert_eq!(actual, 3)
    }

    #[test]
    fn test_counting_group_galaxies() {
        let groups = [1, 1, 2, 2, -1];
        let another_groups = [100, 100, 100, 20, 20, 30, 30, -1, 50, 50];
        assert_eq!(count_group_galaxies(&groups, -1), 4);
        assert_eq!(count_group_galaxies(&another_groups, -1), 9);
    }

}