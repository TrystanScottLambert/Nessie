use extendr_api::prelude::*;
use ndarray::Array2;

#[derive(Debug)]
struct FoFResult {
    pairs: Vec<(usize, usize)>,
    radial_distances: Vec<f64>,
    redshift_differences: Vec<f64>,
    truncated: bool,
}

fn fof_linking(
    positions: &Array2<f64>,
    bgal: &[f64],
    zcut: &[f64],
    nmax: usize,
) -> FoFResult {
    let n = positions.nrows();

    let mut pairs = Vec::with_capacity(nmax.min(n * (n - 1) / 2));
    let mut radial_distances = Vec::with_capacity(nmax.min(n * (n - 1) / 2));
    let mut redshift_differences = Vec::with_capacity(nmax.min(n * (n - 1) / 2));

    for i in 0..n {
        for j in (i + 1)..n {
            let z_thresh = 0.5 * (zcut[i] + zcut[j]);
            let z_diff = (positions[[i, 3]] - positions[[j, 3]]).abs();

            if z_diff > z_thresh {
                continue;
            }

            let b2 = 0.5 * (bgal[i] + bgal[j]);
            let b2 = b2 * b2;

            let r2 = (0..3)
                .map(|k| (positions[[i, k]] - positions[[j, k]]).powi(2))
                .sum::<f64>();

            if r2 <= b2 {
                if pairs.len() >= nmax {
                    return FoFResult {
                        pairs,
                        radial_distances,
                        redshift_differences,
                        truncated: true,
                    };
                }

                pairs.push((i + 1, j + 1)); // 1-based indexing
                radial_distances.push(r2.sqrt());
                redshift_differences.push(z_diff);
            }
        }
    }

    FoFResult {
        pairs,
        radial_distances,
        redshift_differences,
        truncated: false,
    }
}

/// r entry function to the rust fof algorithm
/// @export
#[extendr]
fn fof_links(ra: Vec<f64>, dec: Vec<f64>, redshifts: Vec<f64>, ) -> list {

}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod Nessie;
    fn fof_links;
}