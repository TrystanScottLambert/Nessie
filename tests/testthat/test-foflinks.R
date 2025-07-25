test_that("On sky differences are caught", {
  ras <-  c(121.1, 121.2, 121.3, 181.1, 181.2, 181.1)
  decs <-  c(-23.1, -23.1, -23.1, 68., 68., 180.)
  distances <-  c(20., 20., 20., 20., 20., 20.)
  linking_lengths <-  c(2., 2., 2., 2., 2., 2.)
  #with the above settings b0 = tan(0.2deg)*10 = 0.0349
  #should connect [(0:2), (3:4), (5)]
  b0 = rep(0.0349, length(ras))
  r0 = rep(123., length(ras)) # high r0 to show isolate to on-sky

  result <-  fof_links_aaron(ras, decs, distances, b0, r0)
  expect_equal(result$i, c(1, 1, 2, 4))
  expect_equal(result$j, c(2, 3, 3, 5))
})

test_that("Line of sight differences are caught", {
  ras <-  c(121.1, 121.2, 121.3, 181.1, 181.2, 181.1)
  decs <-  c(-23.1, -23.1, -23.1, 68., 68., 68.)
  distances <-  c(20., 20., 20., 2., 2., 2.1)
  linking_lengths <-  c(2., 2., 2., 2., 2., 2.)
  #with the above settings r0 = 0.1/2 = 0.05 when b0=1 and Delta Dc = 0.1
  #should connect [(0:2), (3:4), (5)]
  b0 = rep(1., length(ras)) # fixing b0 to isolate line-of-sight
  r0 = rep(0.049, length(decs))
  result <- fof_links_aaron(ras, decs, distances, b0, r0)
  expect_equal(result$i, c(1, 1, 2, 4))
  expect_equal(result$j, c(2, 3, 3, 5))
})

