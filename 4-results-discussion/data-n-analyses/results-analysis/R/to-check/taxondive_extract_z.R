taxondive_extract_z <- function(taxondive.res) {
  ## summarizes the output from the vegan taxondive function and extracts the z
  ## test on the Delta + value and the corresponding p-value to a data frame. 
  
  # get the summary of the taxondive.res object (includes a z test on Delta +, which 
  # we want)
  current.summary <- summary(taxondive.res)
  
  # get the names of the columns, which are in the second part of the summary 
  # object.
  current.names <- dimnames(current.summary)[[2]]
  
  # extract the values of z and the corresponding p-values (found in the 
  # last-but-one and the last columns, respectively). Also, get rid of the last 
  # row, which contains the expected values and returns NA. 
  z.Dplus <- z[1:length(taxondive.res[[1]]), length(current.names) - 1]
  p.value <- z[1:length(taxondive.res[[1]]), length(current.names)]
  
  # round the numbers to a more manageable size
  z.Dplus <- round(z.Dplus, 4)
  p.value <- round(p.value, 6)
  
  # make a data frame for the output
  z.test <- cbind(z.Dplus, p.value)

  return(z.test)
}
