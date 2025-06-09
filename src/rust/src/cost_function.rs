use std::collections::{HashMap, HashSet};


pub fn count_bijective_groups(fof_groups: &[i32], mock_groups: &[i32], singleton_id: i32) -> usize {
    let mut group1_to_indices: HashMap<i32, Vec<usize>> = HashMap::new();
    let mut group2_to_indices: HashMap<i32, Vec<usize>> = HashMap::new();

    for (i, (&g1, &g2)) in fof_groups.iter().zip(mock_groups.iter()).enumerate() {
        if g1 != singleton_id {
            group1_to_indices.entry(g1).or_default().push(i);
        }
        if g2 != singleton_id {
            group2_to_indices.entry(g2).or_default().push(i);
        }
    }

    let mut matched_mock_groups: HashSet<i32> = HashSet::new();
    let mut bijective_matches = 0;

    for indices_fof in group1_to_indices.values() {
        let mut best_match: Option<i32> = None;

        for &idx in indices_fof {
            let mock_id = mock_groups[idx];
            if mock_id == singleton_id || matched_mock_groups.contains(&mock_id) {
                continue;
            }

            if let Some(indices_mock) = group2_to_indices.get(&mock_id) {
                let set_fof: HashSet<_> = indices_fof.iter().copied().collect();
                let set_mock: HashSet<_> = indices_mock.iter().copied().collect();
                let overlap: HashSet<_> = set_fof.intersection(&set_mock).copied().collect();

                let overlap_len = overlap.len() as f64;
                let fof_len = set_fof.len() as f64;
                let mock_len = set_mock.len() as f64;

                if overlap_len / fof_len > 0.5 && overlap_len / mock_len > 0.5 {
                    best_match = Some(mock_id);
                    break;  // Stop after first bijective match
                }
            }
        }

        if let Some(mock_id) = best_match {
            bijective_matches += 1;
            matched_mock_groups.insert(mock_id); // Prevent duplicate use
        }
    }

    bijective_matches
}

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
pub fn calculate_e_metrics(fof_groups: &[i32], mock_groups: &[i32], singleton_id: i32) -> QuantityResult {
    let n_bij = count_bijective_groups(fof_groups, mock_groups, singleton_id) as f64;
    let n_fof = n_bij / count_groups(fof_groups, -1) as f64;
    let n_mock = n_bij / count_groups(mock_groups, -1) as f64;
    QuantityResult { fof_metric: n_fof, mock_metric: n_mock, total_metric: n_fof * n_mock }
}


#[derive(Debug, Clone)]
pub struct GroupCalcs {
    pub purity: f64,
    pub length: usize
}

pub fn calculate_purity_per_group(
    group_cat_1: &[i32],
    group_cat_2: &[i32],
    singleton_id: i32,
) -> HashMap<i32, GroupCalcs> {
    let mut group1_to_indices: HashMap<i32, Vec<usize>> = HashMap::new();
    let mut group2_to_indices: HashMap<i32, Vec<usize>> = HashMap::new();

    for (i, (&g1, &g2)) in group_cat_1.iter().zip(group_cat_2.iter()).enumerate() {
        if g1 != singleton_id {
            group1_to_indices.entry(g1).or_default().push(i);
        }
        if g2 != singleton_id {
            group2_to_indices.entry(g2).or_default().push(i);
        }
    }

    let mut results: HashMap<i32, GroupCalcs> = HashMap::new();

    for (&group_id, indices_in_fof) in &group1_to_indices {
        let mut best_score: f64 = 0.0;

        let mut seen_mock_groups = HashSet::new();
        for &idx in indices_in_fof {
            let g2 = group_cat_2[idx];
            if g2 != singleton_id {
                seen_mock_groups.insert(g2);
            }
        }

        for &mock_group in &seen_mock_groups {
            let indices_in_mock = &group2_to_indices[&mock_group];

            let overlap: HashSet<_> = indices_in_fof
                .iter()
                .copied()
                .collect::<HashSet<_>>()
                .intersection(&indices_in_mock.iter().copied().collect())
                .copied()
                .collect();

            let n_common = overlap.len() as f64;
            let n_fof = indices_in_fof.len() as f64;
            let n_mock = indices_in_mock.len() as f64;

            let score = (n_common / n_fof) * (n_common / n_mock);
            best_score = best_score.max(score);
        }

        results.insert(
            group_id,
            GroupCalcs {
                purity: best_score,
                length: indices_in_fof.len(),
            },
        );
    }

    results
}

pub fn calculate_q_score(group_cat_1: &[i32], group_cat_2: &[i32], singleton_id: i32) -> f64 {
    let numerator = calculate_purity_per_group(group_cat_1, group_cat_2, singleton_id);
    let denominator = count_group_galaxies(group_cat_1, singleton_id) as f64;
    let sum: f64 = numerator.values().map(|calcs| (calcs.length as f64) * calcs.purity).sum();
    sum / denominator
}

pub fn calculate_q_metrics(fof_groups: &[i32], mock_groups: &[i32], singleton_id: i32) -> QuantityResult {
    let q_fof = calculate_q_score(fof_groups, mock_groups, singleton_id);
    let q_mock = calculate_q_score(mock_groups, fof_groups, singleton_id);
    QuantityResult {fof_metric: q_fof, mock_metric: q_mock, total_metric: q_fof * q_mock}
}

pub struct STotalMetrics {
    pub e_metrics: QuantityResult,
    pub q_metrics: QuantityResult,
    pub s_total: f64
}

pub fn calculate_s_total(fof_groups: &[i32], mock_groups: &[i32], singleton_id: i32) -> STotalMetrics {
    let e_metrics = calculate_e_metrics(fof_groups, mock_groups, singleton_id);
    let q_metrics = calculate_q_metrics(fof_groups, mock_groups, singleton_id);
    let s_total = e_metrics.total_metric * q_metrics.total_metric;
    STotalMetrics { e_metrics, q_metrics, s_total}
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
    fn test_bij_perfect_score() {
        // ideal case where everything is perfectly matched.
        let fof_groups = [1, 1, 1, 2, 2, 2];
        let mock_groups = [1, 1, 1, 2, 2, 2];
        let actual = count_bijective_groups(&fof_groups, &mock_groups, -1);
        let expected = 2;
        assert_eq!(actual, expected);

        // shouldn't matter if labels are switched around either
        let fof_groups = [1, 1, 1, 2, 2, 2];
        let mock_groups = [2, 2, 2, 1, 1, 1];
        let actual = count_bijective_groups(&fof_groups, &mock_groups, -1);
        let expected = 2;
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_bij_nothing() {
        // the case where there are no matches.
        let fof_groups = [1, 1, 1, 2, 2, 2];
        let mock_groups = [1, 2, 3, 3, 5, 6];
        let actual = count_bijective_groups(&fof_groups, &mock_groups, -1);
        let expected = 0;
        assert_eq!(actual, expected);

        // case with isolated galaxies in the fof groups
        let fof_groups = [1, 1, 1, -1, -1, -1];
        let mock_groups = [99, 98, 97, 5, 5, 5];
        let actual = count_bijective_groups(&fof_groups, &mock_groups, -1);
        let expected = 0;
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_bij_complicated() {
        let fof_groups = [1, 1, 1, 2, 2, -1, -1];
        let mock_groups = [88, 88, 90, 90, 90, 100, 100];
        let actual = count_bijective_groups(&fof_groups, &mock_groups, -1);
        let expected = 2;
        assert_eq!(actual, expected);

        // Inverse it. Shouldn't matter.
        let actual = count_bijective_groups(&mock_groups, &fof_groups, -1);
        let expected = 2;
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_bij_larger_fof_group() {
        // testing the case where the fof group is larger than the mock group.s
        let fof_groups = [-1, 1, 1, 1, 1];
        let mock_groups = [-1, 100, 100, -1, -1];
        let actual = count_bijective_groups(&fof_groups, &mock_groups, -1);
        assert_eq!(actual, 0);
    }

    #[test]
    fn test_e_metric() {
        let fof_groups = [1, 1, 1, 2, 2, 3, 3, 3, -1, -1];
        let mock_groups  = [80, 80, 90, 90, 90, 180, 180, 180, 180, -1];
        let actual = calculate_e_metrics(&fof_groups, &mock_groups, -1);
        assert_eq!(actual.total_metric, 1.);

        // only one group ((1/1) * (1/3) = 1/3)
        let fof_groups = [1, 1, -1, -1, -1, -1, -1, -1];
        let mock_groups = [80, 80, 80, 90, 90, 180, 180, -1];
        let actual = calculate_e_metrics(&fof_groups, &mock_groups, -1);
        assert_eq!(actual.total_metric, 1./3.)
    }

    #[test]
    fn test_counting_group_galaxies() {
        let groups = [1, 1, 2, 2, -1];
        let another_groups = [100, 100, 100, 20, 20, 30, 30, -1, 50, 50];
        assert_eq!(count_group_galaxies(&groups, -1), 4);
        assert_eq!(count_group_galaxies(&another_groups, -1), 9);
    }

    #[test]
    fn test_purity_score_from_paper() {
        // using the example in Robotham+2011
        let fof_group = [-1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1];
        let mock_group = [1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2];
        let actual = calculate_purity_per_group(&fof_group, &mock_group, -1);
        assert_eq!(actual.len(), 1);
        let expected = 0.2666667;
        assert!((actual[&1].purity - expected).abs()  < 1e-5);
        assert_eq!(actual[&1].length, 5);
    }

    #[test]
    fn test_purity_score_complicated() {
        let fof_groups = [1, 1, 1, 2, 2, 2, 2, 3, 3, -1, -1];
        let mock_groups = [100, 100, 100, 200, 200, -1, -1, -1, -1, -1, -1];
        let actual = calculate_purity_per_group(&fof_groups, &mock_groups, -1);
        assert_eq!(actual.len(), 3);

        assert_eq!(actual[&1].length, 3);
        assert_eq!(actual[&2].length, 4);
        assert_eq!(actual[&3].length, 2);

        assert_eq!(actual[&1].purity, 1.);
        assert_eq!(actual[&2].purity, 0.5);
        assert_eq!(actual[&3].purity, 0.);


        // inverting should give different answers
        let actual = calculate_purity_per_group(&mock_groups, &fof_groups, -1);
        
        assert_eq!(actual.len(), 2);
        assert_eq!(actual[&100].length, 3);
        assert_eq!(actual[&200].length, 2);

        assert_eq!(actual[&100].purity, 1.);
        assert_eq!(actual[&200].purity, 0.5);

    }

    #[test]
    fn test_calculate_q_score() {
        let fof_groups = [1, 1, 1, 2, 2, 2, 2, 3, 3, -1, -1];
        let mock_groups = [100, 100, 100, 200, 200, -1, -1, -1, -1, -1, -1];

        let q_fof = ((3. * 1.) + (4. * 0.5) + (0.))/9.;
        let q_mock = ((3. * 1.) + (2. * 0.5))/5.;
        let q_total = q_fof * q_mock;

        let actual_q = calculate_q_metrics(&fof_groups, &mock_groups, -1);
        assert!((actual_q.fof_metric - q_fof).abs() < 1e-6);
        assert!((actual_q.mock_metric - q_mock).abs() < 1e-6);
        assert!((actual_q.total_metric - q_total).abs() < 1e-6);
    }

    #[test]
    fn test_calculate_s_total() {
        let fof_groups = [1, 1, 1, 2, 2, 2, 2, 3, 3, -1, -1];
        let mock_groups = [100, 100, 100, 200, 200, -1, -1, -1, -1, -1, -1];

        let n_bijective = 1.;
        let n_fof = 3.;
        let n_mock = 2.;
        let e_fof = n_bijective/n_fof;
        let e_mock = n_bijective/n_mock;
        let e_total = e_fof * e_mock;

        let q_fof = ((3. * 1.) + (4. * 0.5) + (0.))/9.;
        let q_mock = ((3. * 1.) + (2. * 0.5))/5.;
        let q_total = q_fof * q_mock;

        let s_total = q_total * e_total;
        let actual_s = calculate_s_total(&fof_groups, &mock_groups, -1);
        assert!((actual_s.q_metrics.fof_metric - q_fof).abs() < 1e-6);
        assert!((actual_s.q_metrics.mock_metric - q_mock).abs() < 1e-6);
        assert!((actual_s.q_metrics.total_metric - q_total).abs() < 1e-6);
        
        assert!((actual_s.e_metrics.mock_metric - e_mock).abs() < 1e-6);
        assert!((actual_s.e_metrics.fof_metric - e_fof).abs() < 1e-6);
        assert!((actual_s.e_metrics.total_metric - e_total).abs() < 1e-6);
        assert!((actual_s.s_total - s_total).abs() < 1e-6);
    }

}