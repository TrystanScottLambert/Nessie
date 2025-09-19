library(bit64) # needed for true 64-bit integers

test_that("remap_ids works for strings", {
  bad <- c("2mass1", "2mass2", -1, "2mass1", "2mass3")
  good <- c(1, 2, -1, 1, 3)
  ans <- remap_ids(bad)
  expect_equal(ans, good)
})

test_that("remap_ids works for floats", {
  bad <- c(1.1, 2.2, -1, 1.1, 3.3)
  good <- c(1, 2, -1, 1, 3)
  ans <- remap_ids(bad)
  expect_equal(ans, good)
})

test_that("remap_ids works for shark ids", {
  shark_ids <- as.integer64(c(
    21826700000225, 21826700000235, -1,
    21826700000225, 21826700000525
  ))
  good <- c(1, 2, -1, 1, 3)
  ans <- remap_ids(shark_ids)
  expect_equal(ans, good)
})
