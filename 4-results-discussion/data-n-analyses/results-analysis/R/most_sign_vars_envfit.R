sign_vars_pos_envfit <- function(envfit.obj) {
  ## helper function to extract the most significant (p < 0.05) variables
  ## and their position (based on p and r2).
  
  # extract the most significant variables (p < 0.05) from the fit
  envfit.sign.vars <- extract_envfit_scores(envfit.obj, p = 0.05, r2 = TRUE)
  
  # sort the resulting data frame by p-value and descending r-squared
  envfit.sign.vars <- arrange(envfit.sign.vars, pvals, desc(r2))
  
  # count the number of times a variable occurs at a particular position:   
  # only get vars + row names (variable & position) from the data frame
  envfit.sign.vars$pos <- row.names(envfit.sign.vars)
  sign.vars.pos <- subset(envfit.sign.vars, select = c(vars, pos))
  
  # covert pos to numeric for easier sorting
  sign.vars.pos$pos <- as.numeric(sign.vars.pos$pos)
  
  return(sign.vars.pos)  
} 


sign_vars_freq <- function(env.vars.count, target.freq) {
  ## completely useless helper function, to be applied only here, for examining
  ## the frequency of occurence of environemental variables across a list of 
  ## variables and their positions obtained from the respective envfits. 
  ## Extracts the names of the final variables - to be plotted as surfaces or
  ## vectors over an mds.
  ## Input is a list of data frames (each - columns vars and pos).
  
  # convert list of data frames to 1 big data frame
  sign.count.df <- do.call("rbind", env.vars.count)
  
  # 1. count the total number of times a variable appears as significant across the 
  # imputed datasets
  sign.count.all <- arrange(count(sign.count.df, vars = "vars"), desc(freq))
  
  # 2. count the occurrences of each variable in each position
  sign.count.pos <- count(sign.count.df, vars = c("vars", "pos"))
  sign.count.pos <- arrange(sign.count.pos, pos, desc(freq))
  
  # examine the distribution of the variables and select the most significant 
  # (the ones that appear the most often, i.e. with a frequency >= target frequency)
  # at the topmost positions after the envfit procedure)
  sign.vars <- unique(subset(sign.count.pos, freq >= target.freq)$vars)
  
  # drop the unused factor levels & reorder the factor (order is important here = 
  # decreasing significance)
  sign.vars <- droplevels(sign.vars)
  sign.vars <- reorder.factor(sign.vars, new.order = unique(sign.vars))
  
  # return the variables as well as all counts in a list, just in case
  sign.vars.list <- list("sign.vars" = sign.vars, 
                         "sign.freq" = sign.count.all, 
                         "sign.pos.freq" = sign.count.pos) 
  return(sign.vars.list)
}
