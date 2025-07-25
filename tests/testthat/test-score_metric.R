test_that("score is the same as the old code", {
  library(arrow)

  # Setting up cosmology
  cosmo <- FlatCosmology$new(h = 1, omega_matter = 0.3)
  g09_randoms <- as.data.frame(read.csv(
    "~/Desktop/GAMA_paper_plotter/gama_g09_randoms.txt"
  ))
  g09_rho_mean <- Nessie::create_density_function(
    g09_randoms$z,
    length(g09_randoms$z) / 400,
    0.001453924,
    cosmo
  )

  # Setting up Redshift catalogues
  g09_lightcone <- as.data.frame(arrow::read_parquet(
    "~/Desktop/GAMA_paper_plotter/mocks/gama_gals_for_R.parquet"
  ))
  g09_lightcone <- g09_lightcone[
    g09_lightcone['lightcone_gamafield'] == '0g09',
  ]
  g09_lightcone <- g09_lightcone[g09_lightcone["zobs"] < 0.5, ]

  red_cat <- RedshiftCatalog$new(
    g09_lightcone$ra,
    g09_lightcone$dec,
    g09_lightcone$zobs,
    g09_rho_mean,
    cosmo
  )
  red_cat$set_completeness()

  #Annoyingly we have to change the -1 thing to the old GroupID method to to do this test.
  neg1_idx <- which(g09_lightcone$GroupID == -1)
  maxgnum <- max(g09_lightcone$GroupID[g09_lightcone$GroupID != -1])
  g09_lightcone$GroupID[neg1_idx] <- seq(
    from = maxgnum + 1,
    length.out = length(neg1_idx)
  )
  red_cat$mock_group_ids <- g09_lightcone$GroupID

  b0 <- 0.05
  r0 <- 18
  ### old method from GAMA group finder
  grefs <- red_cat$get_raw_groups(b0, r0)
  maxgnum <- max(grefs[, 2])
  singles <- which(!seq_along(g09_lightcone$ra) %in% grefs[, 1])
  singletons <- cbind(singles, (maxgnum + 1):(maxgnum + length(singles)))
  colnames(singletons) <- colnames(grefs)
  fullgroup <- rbind(grefs, singletons)
  fullgroup <- cbind(fullgroup$galaxy_id, fullgroup$group_id)

  start.now <- Sys.time()
  mock_grefs <- cbind(seq_along(g09_lightcone$ra), g09_lightcone$GroupID)
  mock_t <- .bijcheck(mock_grefs, fullgroup, groupcut = 2)
  fof_t <- .bijcheck(fullgroup, mock_grefs, groupcut = 2)
  score <- (fof_t$summary['bij'] * mock_t$summary['bij']) *
    (fof_t$summary['int'] * mock_t$summary['int'])

  mock_t <- .bijcheck(mock_grefs, fullgroup, groupcut = 3)
  fof_t <- .bijcheck(fullgroup, mock_grefs, groupcut = 3)
  score_3 <- (fof_t$summary['bij'] * mock_t$summary['bij']) *
    (fof_t$summary['int'] * mock_t$summary['int'])
  end.now <- Sys.time()
  print(end.now - start.now)
  ###

  ### new method
  ids <- g09_lightcone$GroupID
  counts <- table(ids)
  singleton_ids <- names(counts[counts == 1])
  red_cat$mock_group_ids <- ifelse(ids %in% singleton_ids, -1, ids)
  red_cat$run_fof(b0, r0)
  start.now <- Sys.time()
  score_me <- red_cat$compare_to_mock()
  end.now <- Sys.time()
  print(end.now - start.now)

  expect_equal(score_me, as.numeric(score), tolerance = 1e-2)
  expect_equal(
    red_cat$compare_to_mock(3),
    as.numeric(score_3),
    tolerance = 1e-2
  )
})
