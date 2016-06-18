alpha_diversity <- function(zoo.data) {
  ## uses vegan package functions to calculate species richness S, Shannon H', 
  ## Pielou J', and effective Shannon diversity (the Hill number), and returns 
  ## them in a data frame.
  ## Arguments: zoo.data - abundance data; samples x species (vegan format)

  library(vegan)
  
  # make a subset of the numeric columns (the species); exclude the factors
  # (station names, replicates, etc.)
  num.zoo <- zoo.data[sapply(zoo.data, is.numeric)]
  
  # number of species per sample/replicate
  s <- specnumber(num.zoo)
  
  # Shannon-Wiener index - only on the numeric columns (species)! -> quantifies 
  # the probability that any two species drawn from the community are different - 
  # measure of entropy, NOT diversity! - doesn't behave linearly! Actually measures
  # the level of disorder in the system (more disorder = more diversity).
  h <- diversity(num.zoo)
  
  # Pielou evenness - level of dominance in the community. Close to 1 - abundance 
  # more evenly distrbuted among the species; close to 0 - one or several species 
  # largely dominate in abundance over the others.
  j <- h/log(s) 
  
  # effective diversity - to allow a meaningful comparison of diversity, not of 
  # entropy (=Shannon, Simpson, Gini-Simpson) -> Hill numbers. 
  # Effective numbers impose a linear transformation on Shannon diversity - 
  # represent the number of equally abundant species necessary to produce the observed
  # value of diversity. Range from 1 to S, where a value of S would indicate all 
  # species are present and in equal abundances.
  # The conversion of Shannon diversity to effective numbers is exp(H)
  h.effective <- exp(h)
  
  # make a data frame of the different diversity indices and the factors describing
  # the samples (station name, date, replicate..)
  diversity.table <- cbind(zoo.data[sapply(zoo.data, is.factor)], 
                           s, 
                           h, 
                           j, 
                           h.effective)
  
  return(diversity.table)
}