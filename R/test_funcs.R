# functions used to test against the previous version

#' @export
bijcheck = function(grefs1, grefs2, groupcut = 2) {
  tab1 = table(grefs1[, 2])
  tab2 = table(grefs2[, 2])
  tab1 = cbind(as.numeric(names(tab1)), as.numeric(tab1))
  tab2 = cbind(as.numeric(names(tab2)), as.numeric(tab2))
  tab1 = tab1[tab1[, 2] >= groupcut, ]
  #tab2 = tab2[tab2[, 2]>1, ]
  grefs1 = grefs1[grefs1[, 2] %in% tab1[, 1], ]
  grefs2 = grefs2[grefs2[, 2] %in% tab2[, 1], ]
  grouplist = sort(unique(grefs1[, 2]))
  bijout = {}
  for (i in seq_along(grouplist)) {
    tempgrps = grefs2[grefs2[, 1] %in% grefs1[grefs1[, 2] == grouplist[i], 1], 2]
    temptab = table(tempgrps)
    temptab = cbind(as.numeric(names(temptab)), as.numeric(temptab))
    frac1 = temptab[, 2] / tab1[tab1[, 1] == grouplist[i], 2]
    frac2 = temptab[, 2] / tab2[tab2[, 1] %in% temptab[, 1], 2]
    add = tab1[tab1[, 1] == grouplist[i], 2] - length(tempgrps)
    if (length(add) > 0) {
      if (add < 0) {
        print(add)
        print(i)
      }
      temptab = rbind(temptab, cbind(rep(NA, add), rep(1, add)))
      frac1 = c(frac1, rep(1, add) / tab1[tab1[, 1] == grouplist[i], 2])
      frac2 = c(frac2, rep(1, add))
    }
    comb = frac1 * frac2
    frag = length(comb)
    bestloc = which.max(comb)
    #Assuming both the pregenerated catalogue and the fof catalogue are sorted by group ID (true generally) matchrow
    #will give the row ref of the matching group, saves doing this externally
    matchrow = which(tab2[, 1] == temptab[bestloc, 1])
    if (length(matchrow) == 0) {
      matchrow = NA
    }
    if (is.na(temptab[bestloc, 1]) == FALSE) {
      bijout = rbind(bijout, c(grouplist[i], temptab[bestloc, 1], i, matchrow,
                               tab1[tab1[, 1] == grouplist[i], 2], tab2[tab2[, 1] == temptab[bestloc, 1], 2],
                               frag, frac1[bestloc], frac2[bestloc], comb[bestloc]))
    } else {
      bijout = rbind(bijout, c(grouplist[i], NA, i, NA, tab1[tab1[, 1] == grouplist[i], 2], 1, frag,
                               frac1[bestloc], frac2[bestloc], comb[bestloc]))
    }
  }
  colnames(bijout) = c("GID1", "GID2", "RID1", "RID2", "N1", "N2", "Frag", "Q1", "Q2", "TotQ")
  #Below I generate some summary statistics
  G1bij_num = length(which(bijout[, "Q1"] > 0.5 &  bijout[, "Q2"] > 0.5))
  G1bij_den = length(bijout[, "N1"])
  G1bij = G1bij_num / G1bij_den
  G1perfect = length(which(bijout[, "Q1"] == 1 &  bijout[, "Q2"] == 1)) / length(bijout[, "N1"])
  G1int_num = sum(bijout[, "Q1"] * bijout[, "N1"])
  G1int_den = sum(bijout[, "N1"])
  G1int = G1int_num / G1int_den
  G1eff = G1bij * G1int
  return(list(groups = bijout, summary = c(bij_num = G1bij_num, bij_den = G1bij_den, bij = G1bij, per = G1perfect,
                                           int_num = G1int_num, int_den = G1int_den, int = G1int, eff = G1eff)))
}

#' @export
bijcheck_simple = function(grefs1, grefs2, groupcut = 2) {
  # Create frequency tables
  tab1 = table(grefs1[, 2])
  tab2 = table(grefs2[, 2])
  tab1 = cbind(as.numeric(names(tab1)), as.numeric(tab1))
  tab2 = cbind(as.numeric(names(tab2)), as.numeric(tab2))

  # Filter groups by size threshold
  tab1 = tab1[tab1[, 2] >= groupcut, ]
  grefs1 = grefs1[grefs1[, 2] %in% tab1[, 1], ]
  grefs2 = grefs2[grefs2[, 2] %in% tab2[, 1], ]

  grouplist = sort(unique(grefs1[, 2]))

  # Initialize vectors to store Q1 values and N1 values
  Q1_values = numeric(length(grouplist))
  N1_values = numeric(length(grouplist))
  Q2_values = numeric(length(grouplist))

  for (i in seq_along(grouplist)) {
    # Find overlapping items between current group and grefs2
    tempgrps = grefs2[grefs2[, 1] %in% grefs1[grefs1[, 2] == grouplist[i], 1], 2]

    if (length(tempgrps) > 0) {
      temptab = table(tempgrps)
      temptab = cbind(as.numeric(names(temptab)), as.numeric(temptab))

      # Calculate fractions
      N1_current = tab1[tab1[, 1] == grouplist[i], 2]
      frac1 = temptab[, 2] / N1_current
      frac2 = temptab[, 2] / tab2[tab2[, 1] %in% temptab[, 1], 2]

      # Find best match
      comb = frac1 * frac2
      bestloc = which.max(comb)

      Q1_values[i] = frac1[bestloc]
      Q2_values[i] = frac2[bestloc]
    } else {
      Q1_values[i] = 0
      Q2_values[i] = 0
    }

    N1_values[i] = tab1[tab1[, 1] == grouplist[i], 2]
  }

  # Calculate the four required values
  G1bij_num = sum(Q1_values > 0.5 & Q2_values > 0.5)
  G1bij_den = length(N1_values)
  G1int_num = sum(Q1_values * N1_values)
  G1int_den = sum(N1_values)

  return(list(
    G1bij_num = G1bij_num,
    G1bij_den = G1bij_den,
    G1int_num = G1int_num,
    G1int_den = G1int_den
  ))
}

#' @export
bijcheck_simple_minus1 = function(grefs1, grefs2, groupcut = 2) {
  # Create frequency tables (this will automatically exclude -1 values)
  tab1 = table(grefs1[grefs1 != -1])
  tab2 = table(grefs2[grefs2 != -1])
  tab1 = cbind(as.numeric(names(tab1)), as.numeric(tab1))
  tab2 = cbind(as.numeric(names(tab2)), as.numeric(tab2))

  # Filter groups by size threshold
  tab1 = tab1[tab1[, 2] >= groupcut, ]

  # Get indices of galaxies in valid groups
  valid_indices1 = which(grefs1 %in% tab1[, 1])
  valid_indices2 = which(grefs2 %in% tab2[, 1])

  grouplist = sort(unique(grefs1[valid_indices1]))

  # Initialize vectors to store Q1 values and N1 values
  Q1_values = numeric(length(grouplist))
  N1_values = numeric(length(grouplist))
  Q2_values = numeric(length(grouplist))

  for (i in seq_along(grouplist)) {
    # Find galaxies in current group from grefs1
    group_galaxies = which(grefs1 == grouplist[i])
    # Find what groups those same galaxies belong to in grefs2
    overlap_groups = grefs2[group_galaxies]
    # Remove isolated galaxies from overlaps for table calculation
    overlap_groups_valid = overlap_groups[overlap_groups != -1]

    if (length(overlap_groups_valid) > 0) {
      temptab = table(overlap_groups_valid)
      temptab = cbind(as.numeric(names(temptab)), as.numeric(temptab))

      # Calculate fractions
      N1_current = tab1[tab1[, 1] == grouplist[i], 2]
      frac1 = temptab[, 2] / N1_current
      frac2 = temptab[, 2] / tab2[tab2[, 1] %in% temptab[, 1], 2]

      # Add padding for galaxies that don't have overlaps (are -1 in grefs2)
      num_isolated = sum(overlap_groups == -1)
      if (num_isolated > 0) {
        # Add entries for isolated galaxies
        frac1 = c(frac1, rep(1, num_isolated) / N1_current)
        frac2 = c(frac2, rep(1, num_isolated))
      }

      # Find best match
      comb = frac1 * frac2
      bestloc = which.max(comb)

      Q1_values[i] = frac1[bestloc]
      Q2_values[i] = frac2[bestloc]
    } else {
      # All galaxies in this group are isolated in grefs2
      N1_current = tab1[tab1[, 1] == grouplist[i], 2]
      Q1_values[i] = 1 / N1_current  # Each isolated galaxy contributes 1/N1
      Q2_values[i] = 1  # Each isolated galaxy is 100% of its own "group"
    }

    N1_values[i] = tab1[tab1[, 1] == grouplist[i], 2]
  }

  # Calculate the four required values
  G1bij_num = sum(Q1_values > 0.5 & Q2_values > 0.5)
  G1bij_den = length(N1_values)
  G1int_num = sum(Q1_values * N1_values)
  G1int_den = sum(N1_values)



  return(list(
    G1bij_num = G1bij_num,
    G1bij_den = G1bij_den,
    G1int_num = G1int_num,
    G1int_den = G1int_den
  ))
}

