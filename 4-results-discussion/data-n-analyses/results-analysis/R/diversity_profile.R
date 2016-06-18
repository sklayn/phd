diversity_profiles <- function(community, q = 50) {
  # calculates diversity profiles of a community: effective diversities for a sequence
  # of values of q - the parameter in the general entropy equation that determines 
  # the sensitivity to rare species. High q = more emphasis on abundant species; 
  # low q - more emphasis on rare species; q = 0 - both kinds treated exactly the
  # same (= species richness). 

  # drop the all-0 columns (taxa present in the taxa list, but absent from any of the 
  # current samples)
  community <- community[colSums(community) > 0]
  
  # get the values of q
  q.seq <- seq(length = q, from = 0, by = 0.11)

  # calculate the proportions of each species per sample
  taxa.proportions <- apply(community, 1, function(x) x/sum(x))
  
  # raise only proportions > 0 by power q (otherwise for q = 0 the entropy equation
  # does not correspond to the species richness of the sample); then apply the 
  # general entropy equation for each value of q 

  entropies <- lapply(q.seq, function (q) 
                                colSums(apply(taxa.proportions, 2, 
                                              function(x) ifelse(x > 0, x^q, 0)))^(1/(1-q)))
  
  # clean up the output (namely, get rid of the horrible column names)
  diversity.profiles <- as.data.frame(entropies)
  names(diversity.profiles) <- q.seq
  diversity.profiles <- t(diversity.profiles)
  
  return(diversity.profiles)
}
