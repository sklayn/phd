weighted_div_profiles <- function(community.data, tax.data, q = 50) {
  ## calculates diversity profiles for a community, weighted by a measure of 
  ## similarity between the taxa present (can be phylogenetic distance, a 
  ## functional measure of similarity, etc.). 
  ##
  ## Based on: Leinster & Cobbold, 2012. Measuring diversity: the importance of 
  ## species similarity.
  ##
  ## Arguments: community.data - zoobenthic abundance data, samples x species
  ##            tax.data - taxonomic data (or other similarity measure between
  ##                the species). First column - species name, other columns -
  ##                progressively higher taxonomic categories.
  ##            q - number of values for the sensitivity parameter q, for which
  ##                diversities of order q will be calculated to form a profile.  
  
  library(vegan)
  
  # select only the species present in the current dataset (remove all-0 species)
  community <- community.data[colSums(community.data) > 0]
  
  # get the values of q for which diversities will be calculated
  q.seq <- seq(length = q, from = 0, by = 0.11)
  
  # calculate the proportions of each species per sample; output is species x 
  # samples matrix!
  taxa.proportions <- apply(community, 1, function(x) x/sum(x)) 
  
  # subset the provided taxonomic classification table (or other distance mesure) 
  # by matching species names with community table to only work with the species 
  # present in the current community  
  current.zoo.taxa <- subset(tax.data, row.names(tax.data) %in% names(community))
  
  # compute taxonomic distances using variable step - scaling according to reduction
  # in number of classes when going up in the tax. tree (-> variable - because if 
  # almost all genera have only 1 species, it makes no big difference if 2 individuals 
  # belong to a different species or a different genus).
  taxa.distance <- taxa2dist(current.zoo.taxa, varstep = TRUE)
  
  # convert distance object to matrix, then scale so that distances vary between
  # 0 and 1, where 0 means no difference between taxa, and the dissimilarities 
  # progressively increase (vegan function taxa2dist gives us the opposite).
  taxa.dist.matrix <- as.matrix(taxa.distance)
  taxa.dist.matrix <- 1 - (taxa.dist.matrix / 100)
  
  # calculate the ordinariness of each species within the community: each value is
  # the expected similarity between an individual of the ith species and an individual
  # chosen at random.
  taxa.ordinariness <- taxa.dist.matrix %*% taxa.proportions
  
  weighted_div <- function(prop.mat, ord.mat, q) {
    ## calculates the diversity of order q for a community matrix containing 
    ## species proportions in the samples, weighted by a measure of the 
    ## ordinariness of each species in the community. 
    
    # only calculate for species actually present in the community
    out.mat <- ifelse(ord.mat > 0, ord.mat^(q - 1), 0)  
    
    # calculate the average ordinariness of an individual from the community 
    # (= sum(prop * ord), done by columns because samples/communities are each 
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
                                                           taxa.ordinariness, q))
  
  # convert output list to data frame and clean it up a bit (get rid of
  # the horrible column names and transpose so that samples are columns)
  weighted.profiles <- as.data.frame(weighted.profiles)
  names(weighted.profiles) <- q.seq
  weighted.profiles <- t(weighted.profiles)

  return(weighted.profiles) 
}  

