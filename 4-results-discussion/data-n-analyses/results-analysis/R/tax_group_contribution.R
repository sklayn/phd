tax_group_contribution <- function(abnd.data, tax.groups) {
  ## calculates the proportional contribution of each of a set of larger taxonomic
  ## groups to the overall abundance/biomass in the stations' communities.
  ## Arguments: abnd.data - data frame of numeric abundances/biomasses for all
  ##              taxa present in the current dataset (stations x species)
  ##            tax.data - vector of taxonomic group labels of the same length 
  ##              as the number of non-0 (i.e., present) species in the current
  ##              dataset.
  ## Returns a data frame of proportions of each taxonomic group in each station.
  ## For each station, the sum of all taxonomic groups' proportions (= rows) 
  ## should be 1.
  
  library(plyr)
  
  # calculate each species' proportion in each sample 
  prop.abnd <- abnd.data / rowSums(abnd.data)
  
  # transpose the proportions data frame and add the factor column containing each
  # species' taxonomic group
  prop.abnd <- as.data.frame(t(prop.abnd))
  names(prop.abnd) <- 1:ncol(prop.abnd)
  
  tax.abnd <- cbind("tax.group" = tax.groups, prop.abnd)
  
  # aggregate data by taxonomic group (= sum of proportions of all sp in that group)
  summary.tax.abnd <- ddply(tax.abnd, .(tax.group), colwise(sum, .cols = is.numeric))
  
  # transpose to a final data frame of stations x taxonomic groups 
  summary.tax.abnd <- setNames(as.data.frame(t(summary.tax.abnd[, -1])), summary.tax.abnd[, 1])
  
  return(summary.tax.abnd)
}