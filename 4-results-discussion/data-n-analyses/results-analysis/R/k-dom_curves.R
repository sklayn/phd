abc <- function (abnd.data, biomass.data) {
  ## calculates the cumulative abundance and biomass percentages for plotting ABC 
  ## curves. Arguments: abnd.data, biomass.data - data frames of zoobenthic 
  ## abundance and biomass, respectively (species x stations).

  # sort the abundance and biomass data in decreasing order
  abundance <- sort(abnd.data, decreasing = TRUE)
  biomass <- sort(biomass.data, decreasing = TRUE)
  
  # check if the abundance and biomass are imported correctly - they should 
  # contain the same number of elements unless there is an error in the input 
  # data - e.g. a species for which there is only an abundance or only a 
  # biomass.
  if (sum(abundance != 0) != sum(biomass != 0)) {
    stop("The abundance and biomass vectors have different lengths! Check your input data!")
  }
  
  # get rid of the species absent from the sample 
  abundance <- abundance[abundance > 0]
  biomass <- biomass[biomass > 0]
  
  # make a vector to hold the cumulative abundances
  perc.abnd <- NA
  
  # calculate the percent abundance for the first (most abundant) species
  perc.abnd[1] <- abundance[1] / sum(abundance) * 100
  
  # from the second most abundant species onwards, the percent abundance needs to
  # be added to the previous ones in order to obtain a cumulative %
  for (i in 2:length(abundance)) {
    perc.abnd[i] <- (abundance[i] / sum(abundance)) * 100 + perc.abnd[i-1]
  }
  
  # same procedure for the cumulative biomasses 
  perc.bio <- NA
  perc.bio[1] <- biomass[1] / sum(biomass) * 100
  
  for (i in 2:length(biomass)) {
    perc.bio[i] <- (biomass[i] / sum(biomass)) * 100 + perc.bio[i-1]
  }
  
  # Calculate the difference between the biomass and the abundance curves - needed
  # to quantify the closeness of the curves (= W). 
  
  BiAi <- NA
  for(i in 1:length(perc.abnd)) {
    BiAi[i] <- perc.bio[i] - perc.abnd[i]
  }

  # make a data frame for the ABC results (cumulative abundances and biomasses, 
  # and the difference between them (Bi - Ai))
  abc.table <- as.data.frame(cbind(perc.abnd, perc.bio, BiAi))
  colnames(abc.table) <- c("cumulative.abnd", "cumulative.biomass", "BiAi")

  # calculate the W statistic
  W <- sum(abc.table$BiAi) / (50 * (nrow(abc.table) - 1))
  
  # if there is only one species present, the result will be NaN (= division by
  # 0); correct this by setting W = NA (and you shouldn't be doing this analysis 
  # in that case anyway): 
  if (is.nan(W)) {
    W <- NA
  }
  
  # make a list to hold the results: a table of cumulative abundances and biomasses, 
  # and the W statistic for the replicate/sample
  abc.result <- list()
  abc.result$abc <- abc.table
  abc.result$w <- W
  
  return(abc.result)

}


partial_dominance_curves <- function(abnd.data, biomass.data) {
  ## calculates partial dominance curves from a macrobenthic abundance-biomass
  ## dataset. Arguments: abnd.data, biomass.data - vectors/data frames of abundance
  ## and biomass, respectively.
  
  # get the abundance and biomass data and sort them in decreasing order 
  abundance <- sort(abnd.data, decreasing = TRUE)
  biomass <- sort(biomass.data, decreasing = TRUE)
  
  # check if the abundance and biomass vectors contain the same number of elements
  # (i.e., if there is an error in the input data - a species for which only an
  # abundance or only a biomass value is given)
  if (sum(abundance != 0) != sum(biomass != 0)) {
    stop("The abundance and biomass vectors have different lengths! Check your input data!")
  }
  
  # get rid of double-0 - species absent from the sample, which would lead to a division by 0
  abundance <- abundance[abundance > 0]
  biomass <- biomass[biomass > 0]
  
  # make a vector to hold partial percent abundances 
  perc.abnd <- NA
  # calculate the partial abundances, successively removing each species and calculating 
  # over the rest 
  for (i in 1:length(abundance)) {
    
    if (sum(abundance[i:length(abundance)]) == 0){
      perc.abnd[i] <- NA
    }
    else {
      perc.abnd[i] <- (abundance[i] / sum(abundance[i:length(abundance)])) * 100  
    }    
  }

  # same procedure for the biomass
  perc.bio <- NA
  for (i in 1:length(biomass)) {
    if (sum(biomass[i:length(biomass)]) == 0){
      perc.bio[i] <- NA
    }
    else {
      perc.bio[i] <- (biomass[i] / sum(biomass[i:length(biomass)])) * 100  
    }      
  }
  
  # bind the partial abundance and biomass into a new data frame, without matching
  # the species names to its specific values for partial abundance & biomass.
  abc.partial <- as.data.frame(cbind(perc.abnd, perc.bio))
  colnames(abc.partial) <- c("partial.abnd", "partial.biomass")
  
  return(abc.partial)
}
