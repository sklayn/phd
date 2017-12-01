#### Graphical community analysis functions #### 

#### Diversity profiles #####
diversity_profiles <- function(community.data, q = 50) {
  ## calculates diversity profiles of a community: effective diversities for a 
  ## sequence of values of q - the parameter in the general entropy equation 
  ## that determines the sensitivity to rare species. High q = more emphasis on 
  ## abundant species; low q - more emphasis on rare species; q = 0 - both kinds 
  ## treated exactly the same (= species richness).

  # get the sequence of values of q
  q.seq <- seq(length = q, from = 0, by = 0.11)

  # calculate the proportions of each species per sample
  taxa.proportions <- apply(community.data, 1, function(x) x/sum(x))
  
  # Apply the general entropy equation for each value of q. NB raise only 
  # proportions > 0 by power q (otherwise for q = 0 the entropy equation 
  # does not correspond to the species richness of the sample as it should). 
  entropies <- lapply(q.seq, function (q) 
                                colSums(apply(taxa.proportions, 2, 
                                              function(x) ifelse(x > 0, x^q, 0)))^(1/(1-q)))
  
  # clean up the output (namely, get rid of the horrible column names)
  diversity.profiles <- as.data.frame(entropies)
  names(diversity.profiles) <- q.seq
  
  # convert to data frame and then tibble, which is cleaner and overall more 
  # robust. Make sure to keep the row names (= values of q)
  diversity.profiles <- as.data.frame(t(diversity.profiles))
  diversity.profiles <- tibble::rownames_to_column(diversity.profiles, var = "q")

  return(diversity.profiles)
}


weighted_div_profiles <- function(community.data, dist.data, q = 50) {
  ## calculates diversity profiles for a community, weighted by a measure of 
  ## similarity between the taxa present (can be phylogenetic distance, a 
  ## functional measure of similarity, etc.). 
  ##
  ## Based on: Leinster & Cobbold, 2012. Measuring diversity: the importance of 
  ## species similarity.
  ##
  ## Arguments: community.data - zoobenthic abundance data, samples x species
  ##            dist.data - distance data for the species in the dataset (taxonomic, 
  ##              phylogenetic, or other similarity measure to use for weighting).
  ##            q - number of values for the sensitivity parameter q, for which
  ##              diversities of order q will be calculated to form a profile.  
  
  # get the values of q for which diversities will be calculated
  q.seq <- seq(length = q, from = 0, by = 0.11)
  
  # calculate the proportions of each species per sample; output is species x 
  # samples matrix!
  taxa.proportions <- apply(community.data, 1, function(x) x/sum(x)) 
  
  # convert distance object to matrix, then scale so that distances vary between
  # 0 and 1 (0 = no difference between taxa; dissimilarities progressively 
  # increase (vegan function taxa2dist gives the opposite).
  taxa.dist.matrix <- as.matrix(dist.data)
  taxa.dist.matrix <- 1 - (taxa.dist.matrix / 100)
  
  # calculate the ordinariness of each species within the community: the 
  # expected similarity between an individual of the ith species and an 
  # individual chosen at random.
  taxa.ordinariness <- taxa.dist.matrix %*% taxa.proportions
  
  weighted_div <- function(prop.mat, ord.mat, q) {
    ## calculates the diversity of order q for a community matrix containing 
    ## species proportions in the samples, weighted by a measure of the 
    ## ordinariness of each species in the community. 
    
    # only calculate for species actually present in the community
    out.mat <- ifelse(ord.mat > 0, ord.mat^(q - 1), 0)  
    
    # calculate the average ordinariness of an individual from the community 
    # (= sum(prop * ord) - by columns because samples/communities are each 
    # a column in the input matrix). This quantity is large if most of the 
    # population is concentrated into a few very similar species => the average 
    # ordinariness could be called concentration, and is inversely related to 
    # diversity.
    out.mat <- (colSums(prop.mat * out.mat))^(1/(1 - q))
    
    return(out.mat)
  }
  
  # calculate the weighted diversity profiles for the samples over the sequence 
  # of q values. Varying the parameter q varies the influence on diversity of 
  # ordinary species (those i for which the (Zp)i (= taxa.sim[i,]) is large) 
  # relative to unusual species (those for which it is small).
  weighted.profiles <- lapply(q.seq, function(q) weighted_div(taxa.proportions, 
                                                              taxa.ordinariness, 
                                                              q))
  
  # convert output list to data frame and clean it up a bit (get rid of
  # the horrible column names and transpose so that samples are columns)
  weighted.profiles <- as.data.frame(weighted.profiles)
  names(weighted.profiles) <- q.seq
  weighted.profiles <- as.data.frame(t(weighted.profiles))
  weighted.profiles <- tibble::rownames_to_column(weighted.profiles, var = "q")
  
  return(weighted.profiles) 
}

#### ABC curves ####

