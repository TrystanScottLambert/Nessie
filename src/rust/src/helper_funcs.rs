
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
}
