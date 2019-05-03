normalize_index <- function(index, es.lims) {
  ## computes a weighted index ES normalization, given observed index values and 
  ## ES boundaries for the index (according to van de Bund et al. 2008).
  ## 
  ## index - vector of index EQR values (observed)
  ## es.lims - vector of ES limits for this index's EQRs, going from high to bad.
  
  # make vector of weighting factors for each ES (high to bad). Higher classes have higher weights, and weights are equidistant (by 0.2)
  weights <- seq(0.8, 0, by = -0.2)
  
  # calculate the weighted values with a monstrous ifelse, which is thankfully vectorized
  weighted.vals <- ifelse(index > es.lims[1], weights[1] + 0.2 * ((index - es.lims[1])/(1 - es.lims[1])), # high
                    ifelse(index >= es.lims[2], weights[2] + 0.2 * ((index - es.lims[2])/(es.lims[1] - es.lims[2])), # good
                          ifelse(index >= es.lims[3], weights[3] + 0.2 * ((index - es.lims[3])/(es.lims[2] - es.lims[3])), # moderate
                                 ifelse(index >= es.lims[4], weights[4] + 0.2 * ((index - es.lims[4])/(es.lims[3] - es.lims[4])), # poor
                                        weights[5] + 0.2 * index / es.lims[4])))) # bad
                
  return(weighted.vals)
}


es_qualitative <- function(index, es.lims) {
  ## compares the observed index values with the ES boundaries for the index, 
  ## and returns the qualitative ecoligical state (high to bad) as factor.
  ## 
  ## index - vector of index EQR values (observed)
  ## es.lims - vector of ES limits for this index's EQRs, going from high to bad.
  
  es.qualitative <- ifelse(index > es.lims[1], "high",
  ifelse(index >= es.lims[2], "good",
  ifelse(index >= es.lims[3], "moderate",
  ifelse(index >= es.lims[4], "poor",
  "bad"))))
  
  # turn to factor and order levels 
  es.qualitative <- factor(es.qualitative, levels = c("high", "good", "moderate", "poor", "bad"))
  
  return(es.qualitative)
}


ind_abs_diff <- function(indices, nb.classes = 5) {
  ## calculates the absolute class difference between the supplied metrics 
  ## for the given number of ES classes (according to van de Bund et al. 2008). 
  ##
  ## indices - data frame of observed & normalized index values as columns.
  ## nb.classes - number of ES classes for which to calculate the difference.
  ##
  ## output: matrix of absolute class differences for each pair of input indices.
  
  # assign a class to each index value - n equal classes (1 to n) 
  ind.diff <- apply(indices, 2, function(x) ifelse(x > 1, nb.classes, trunc(x * nb.classes + 1)))
  
  # make pairwise combinations of the columns (= indices). 
  ind.combinations <- combn(ncol(ind.diff), 2)
  # fix the names, to use later for the output df
  colnames(ind.combinations) <- apply(ind.combinations, 2, 
  function(x) paste(names(indices)[x[1]], names(indices)[x[2]], sep = "-")) 
  
  # NB results for each pairwise comparison - returned in rows!
  ind.abs.diff <- plyr::adply(ind.combinations, 2,
  function(x) abs(ind.diff[, x[1]] - ind.diff[, x[2]]))
  
  # remove first column (id, and a useless one at that in this case)
  ind.abs.diff <- ind.abs.diff[-1]
  
  # convert to numeric (everything is a character because of the late id column)
  ind.abs.diff <- apply(ind.abs.diff, 2, as.numeric)
  
  # transpose, so that each comparison is a column
  ind.abs.diff.final <- t(ind.abs.diff)
  
  # fix the column names to make them more informative
  colnames(ind.abs.diff.final) <- colnames(ind.combinations)
  
  return(ind.abs.diff.final)

}
