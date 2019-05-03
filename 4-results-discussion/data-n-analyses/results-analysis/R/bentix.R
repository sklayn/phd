bentix <- function(zoo.data, ecological.groups) { 
  ## function to calculate the biotic index BENTIX
  
  # convert input data to matrix
  zoo.data <- as.matrix(zoo.data)
  
  # assign the species from the input data table to their ecological sensitivity 
  # groups: match the species in the EG table with the species in the input data table.
  eg <- ecological.groups[which(rownames(ecological.groups) %in% colnames(zoo.data)), , drop = F]
  
  # if there are extra species in the input data table (without an EG assignment), set their EG 
  # value to NA and add them to the end of the EG table.
  if(length(which(colnames(zoo.data) %in% rownames(eg))) < length(colnames(zoo.data))) {
    eg.na <- colnames(zoo.data)[-which(colnames(zoo.data) %in% rownames(eg))]
    eg <- rbind(eg, rep(NA, length(eg.na)))
    rownames(eg)[(nrow(eg)-length(eg.na)+1):nrow(eg)] <- colnames(zoo.data)[-which(colnames(zoo.data) %in% rownames(eg))]
  }
  
  # order the EG table alphabetically, to put any newly added species in the correct spot
  eg <- eg[order(rownames(eg)), , drop = F]
  
  # make a new table to hold the EG assignments, which will be used in the BENTIX calculation - 
  # 2 columns for the 2 EGs and an additional column for species without assignment (NAs), and 
  # fill it with 0 only.
  eg.matr <- as.data.frame(matrix(0, nrow = ncol(zoo.data), ncol = 3))
  colnames(eg.matr) <- c("GS", "GT", "NA")
  rownames(eg.matr) <- rownames(eg)
  
  # match the original EG table and the calculation EG table for each species, resulting in a 
  # matrix of 1 and 0: each species is assigned a value of 1 in the column corresponding to its 
  # EG, and 0 in all the other columns.
  options(warn = -1)
  eg.matr[which(eg[, 1] == 1), 1] <- 1
  eg.matr[which(eg[, 1] == 2), 2] <- 1
  eg.matr[which(eg[, 1] == 0), 3] <- 1
  eg.matr[which(is.na(eg[, 1]) == T), 3] <- 1
  options(warn = 0)
  
  # multiply the species abundance input data by the EG calculation matrix to obtain a table of 
  # abundances in each EG (and NA) (columns), for each station/replicate (rows).
  bentix <- as.data.frame(zoo.data %*% as.matrix(eg.matr))
  
  # add in the beginning of this table a column N for the total abundance (N = sum of each row), 
  # and another in the end for BENTIX - for now empty (filled with NAs). 
  bentix <- cbind("N" = apply(zoo.data, 1, sum), bentix, "BENTIX" = rep(NA, nrow(zoo.data)))
  
  # for each station/replicate (= row), calculate the proportion of each EG
  for (i in 1:nrow(zoo.data)) {
    for (j in 2:4) {
      bentix[i, j] <- bentix[i, j] / bentix[i, "N"] * 100
    }
  }
  
  # calculate BENTIX for each station/replicate, rounding the value to 3 digits
  for (i in 1:nrow(zoo.data)) {
    bentix[i, "BENTIX"] <- round((6*bentix[i, "GS"] + 2*bentix[i, "GT"]) / 100, 3)
  }
  
  # rename the EG columns to EG(%); round the proportions to 2 digits
  colnames(bentix)[2:3] <- c("GS(%)", "GT(%)")
  bentix[, 2:4] <- round(bentix[, 2:4], 2)
  
  return(bentix)
}
