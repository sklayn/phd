extract_envfit_scores <- function(envfit.obj, pval = 0.05, r2 = FALSE) {
  # helper function to extract NMDS vector scores from envfit objects according 
  # to the specified p-value (p < 0.05 by default). Optionally also extracts 
  # the r2 value.
  
  # extract the NMDS scores for the fitted environmental parameters 
  vector.scrs <- as.data.frame(scores(envfit.obj, display = "vectors"))
  
  # convert all the variables characterizing the vectors in the envfit object 
  # to a list 
  list.scrs <- as.list(envfit.obj$vectors)
  
  if(r2) {
    # get r2 values 
    vector.scrs$r2 <- list.scrs$r    
  }
  
  # get the p-values and subset the data frame according to them
  vector.scrs$pvals <- list.scrs$pvals
  sign.vectors <- subset(vector.scrs, pvals < pval)
  
  # clean up a little: get rid of the row names - add them as a variable in their
  # own column instead
  sign.vectors <- cbind(vars = rownames(sign.vectors), sign.vectors)
  rownames(sign.vectors) <- 1:nrow(sign.vectors)
  
  return(sign.vectors)
}
