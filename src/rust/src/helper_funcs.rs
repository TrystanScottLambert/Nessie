
// function to mimic the quantile interpolation that R does. 
pub fn quantile_interpolated(sorted: &[f64], quantile: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    let h = (n - 1) as f64 * quantile;
    let i = h.floor() as usize;
    let frac = h - i as f64;

    if i + 1 < n {
        sorted[i] * (1.0 - frac) + sorted[i + 1] * frac
    } else {
        sorted[i]
    }
}


// Harmonic mean
pub fn harmonic_mean(values: Vec<f64>) -> f64 {
    let n = values.iter().len() as f64;
    let summation = values
        .iter()
        .map(|v| 1./v)
        .sum::<f64>();
    n/summation
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn comparing_quantile_to_r() {
        let sorted_array = [1., 2., 3., 4., 5.];
        assert!((quantile_interpolated(&sorted_array, 0.3) - 2.2).abs() < 1e-6);
        assert!((quantile_interpolated(&sorted_array, 0.5) - 3.).abs() < 1e-6);
        assert!((quantile_interpolated(&sorted_array, 0.77) - 4.08).abs() < 1e-6);
        assert!((quantile_interpolated(&sorted_array, 1.) - 5.).abs() < 1e-6);
    }


    #[test]
    fn comparing_to_psych_in_r_harmonic_mean() {
        // Checking that arrays with a zero are zero.
        let values = vec![0., 1., 2., 3.];
        let values_float = vec![0., 0.1, 0., 0.2, 0.3, 0.5];
        let answer = 0.;
        let result = harmonic_mean(values);
        let result_float = harmonic_mean(values_float);
        assert_eq!(result, answer);
        assert_eq!(result_float, answer);

        // Example often used.
        let values = vec![1., 4., 4.];
        let answer = 2.;
        let result = harmonic_mean(values);
        assert_eq!(result, answer);

        // Testing on more realistic floats.
        let values = vec![0.2, 0.3, 0.21, 0.22, 0.31, 0.5, 0.1, 0.12];
        let answers = 0.1941755;
        let result = harmonic_mean(values);
        assert!((answers - result).abs() < 1e-5)


    }
}
