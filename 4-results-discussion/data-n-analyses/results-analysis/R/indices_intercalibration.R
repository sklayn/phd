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

gm_diff <- function(ind.norm, ind.lims) {
  ### helper to find the proportion of cases where the difference of assessments of two indices crosses the good-moderate boundary
  ### Arguments: ind.norm - index assessments (normalized EQRs), in columns.
  ###            ind.lims - limits of the 5 ecological classes for each index, going from high to bad. NB column number and order should match ind.norm! 
  ### Returns: df of the proportion of cases where the disagreement crosses the good-moderate boundary, for each pair of indices 
  ### Dependencies: tidyverse; custom function es_qualitative for turning numeric EQRs into qualitative ES classes
  
  
  ## get the qualitative assessment in 5 classes, according to each index's limits
  ind.qual <- map2(ind.norm, ind.lims, ~es_qualitative(.x, .y))
  
  ## reduce them to two - acceptable (high/good) and unacceptable (moderate to bad)
  level.key <- c(high = "acceptable", good = "acceptable", moderate = "unacceptable", poor = "unacceptable", bad = "unacceptable")
  
  ind.qual <- lapply(ind.qual, function(x) recode(x, !!!level.key))
  
  ## compare the indices two by two. Indices will go in rows here (default method of collapsing a list in ldply - rbind)! 
  ind.comp <- plyr::ldply(ind.qual)
  
  ## construct a matrix of combinations for pairwise comparison - of ROWS.
  ind.combn <- combn(nrow(ind.comp), 2)
  
  ## fix the names, to use later for the output df
  colnames(ind.combn) <- apply(ind.combn, 2, function(x) paste(ind.comp[x[1], 1], ind.comp[x[2], 1], sep = "-")) 
  
  ## compare the indices 2x2 according to the comparison matrix above. Drop the first column (id) from the index df first! 
  gm.diff <- plyr::adply(ind.combn, 2, function(x) abs(ind.comp[-1][x[1], ] - ind.comp[-1][x[2], ]))
  
  ## calculate the proportion of all cases where a disagreement crosses the GM boundary for each pair of indices (1 = disagreement, 0 = agreement)
  gm.diff <- gm.diff %>% 
    mutate(gm_diff = rowSums(select(., -X1))/ncol(gm.diff %>% select(-X1))) %>% 
    select(X1, gm_diff) %>% 
    rename(contrast = X1)
  
  return(gm.diff)
}
