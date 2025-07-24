test_that("calculate_completeness: errors on mismatched input lengths", {
  ra_obs <- c(10, 20, 30)
  dec_obs <- c(10, 20)  # wrong length

  ra_tgt <- c(15, 25)
  dec_tgt <- c(15, 25)

  ra_eval <- c(12, 18)
  dec_eval <- c(12, 18)

  radii <- c(1, 1)

  expect_error(
    calculate_completeness(ra_obs, dec_obs, ra_tgt, dec_tgt, ra_eval, dec_eval, radii),
    "RA and Dec observed arrays are of different lengths"
  )
})

test_that("calculate_completeness: errors on search_radii and eval mismatch", {
  ra_obs <- c(10, 20)
  dec_obs <- c(10, 20)

  ra_tgt <- c(15, 25)
  dec_tgt <- c(15, 25)

  ra_eval <- c(12, 18)
  dec_eval <- c(12, 18)

  radii <- c(1)  # too short

  expect_error(
    calculate_completeness(ra_obs, dec_obs, ra_tgt, dec_tgt, ra_eval, dec_eval, radii),
    "search_radii and eval arrays must have the same dimensions"
  )
})

test_that("calculate_completeness: works on valid input", {
  ra_obs <- c(10, 20)
  dec_obs <- c(10, 20)

  ra_tgt <- c(10, 20, 30)
  dec_tgt <- c(10, 20, 30)

  ra_eval <- c(10, 20)
  dec_eval <- c(10, 20)

  radii <- c(1, 1)

  # Mock the underlying Rust function if not available
  mock_result <- c(1, 0.6667)

  with_mock(
    calc_completeness_rust = function(...) mock_result,
    {
      result <- calculate_completeness(ra_obs, dec_obs, ra_tgt, dec_tgt, ra_eval, dec_eval, radii)
      expect_type(result, "double")
      expect_equal(result, mock_result)
    }
  )
})
