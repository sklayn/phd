#### Graphical community analysis functions #### 

diversity_profiles <- function(community.data, q = 50) {
  ## calculates diversity profiles of a community: 
  ## effective diversities for a sequence of values of q - 
  ## the parameter in the general entropy equation that 
  ## determines the sensitivity to rare species. 
  ## High q = more emphasis on abundant species; low q - 
  ## more emphasis on rare species; q = 0 - both kinds 
  ## treated exactly the same (= species richness).

  # get the sequences of values of q
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


weighted_div_profiles <- function(community.data, tax.data, q = 50) {
  ## calculates diversity profiles for a community, weighted by a measure of 
  ## similarity between the taxa present (can be phylogenetic distance, a 
  ## functional measure of similarity, etc.). 
  ##
  ## Based on: Leinster & Cobbold, 2012. Measuring diversity: the importance of 
  ## species similarity.
  ##
  ## Arguments: community.data - zoobenthic abundance data, samples x species
  ##            tax.data - taxonomic data (or other similarity measure) for 
  ##              the species in the dataset.
  ##            q - number of values for the sensitivity parameter q, for which
  ##              diversities of order q will be calculated to form a profile.  
  
  library(vegan)
  
  # get the values of q for which diversities will be calculated
  q.seq <- seq(length = q, from = 0, by = 0.11)
  
  # calculate the proportions of each species per sample; output is species x 
  # samples matrix!
  taxa.proportions <- apply(community.data, 1, function(x) x/sum(x)) 
  
  # compute taxonomic distances using variable step 
  taxa.distance <- taxa2dist(tax.data, varstep = TRUE)
  
  # convert distance object to matrix, then scale so that distances vary between
  # 0 and 1 (0 = no difference between taxa; dissimilarities progressively 
  # increase (vegan function taxa2dist gives the opposite).
  taxa.dist.matrix <- as.matrix(taxa.distance)
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


plot_div_profiles <- function(div.profiles, stations, one.panel = TRUE) {
  ## plots the diversity profiles by sampling station.
  ## Arguments: div.profiles - matrix of diversity profiles (in wide format). 
  ##              Each profile is a column (rows = q values, columns = 
  ##              samples/stations).
  ##            stations - vector of station names. The order of its levels is 
  ##              assumed to be the correct one, to be used in the plot.
  ##            one.panel - should the profiles for each station be plotted in 
  ##              a single panel (for example when average profiles only), or in
  ##              separate panels by station? Default - single panel.  
  ## Output: ggplot object.
  
  # import the necessary packages
  library(ggplot2)
  library(plyr)
  library(reshape2)  
  
  # reshape data matrix for easier plotting with ggplot
  div.prof.long <- melt(div.profiles)
  names(div.prof.long) <- c("q", "replicate", "value")
  
  # add in the stations to serve as grouping factor
  
  nq <- nrow(div.profiles)                      # number of q values
  nrep <- ncol(div.profiles)/length(stations)   # number of replicates/station
  
  div.prof.long$station <- rep(stations, each = nq * nrep)
  
  # reorder the new factor's levels according to the order of the input factor 
  div.prof.long$station <- reorder.factor(div.prof.long$station,
                                          new.order = stations)
  
  # plot the profiles in one panel (e.g., when one average profile/station)  
  p <- ggplot(data = div.prof.long, aes(x = q, 
                                        y = value, 
                                        group = replicate, 
                                        colour = station)) +
    scale_color_brewer(palette = "Set2") + 
    labs(x = "Sensitivity parameter q", y = "Diversity", colour = "Stations") +
    theme_bw() + 
    theme(legend.text = element_text(size = rel(1))) 
  
  if(!one.panel) {
    # plot the profiles in separate panels by station, if the option is specified
    # at input, and remove the legend 
    p <- p + geom_line() +  
      facet_wrap(~station) +
      theme(legend.position = "none")
  } else {
    # increase the line width to make the lines more visible  
    p <- p + geom_line(size = 0.7)
  }
  
  return(p)
}


plot_div_profiles_w_aver <- function(div.profiles, 
                                     aver.profiles, 
                                     stations, 
                                     col.profiles = c("skyblue", "royalblue")) {
  ## plots the diversity profiles by sampling station, and the average diversity 
  ## profile for each station.
  ## Arguments: div.profiles - matrix of diversity profiles (in wide format). 
  ##              Each profile is a column (rows = q values, columns = 
  ##              samples/stations).
  ##            aver.profiles - matrix of average diversity profiles (in wide 
  ##              format). Each profile is a column (1 per station).
  ##            stations - vector of station names. The order of its levels is 
  ##              assumed to be the correct one, to be used in the plot.
  ##            col.profiles - char vector of colors for the profiles. The first 
  ##              item is used for the full profiles, the second - for the average
  ##              profile.
  ## Output: ggplot object.
  
  # import the necessary libraries
  library(ggplot2)
  library(plyr)
  library(reshape2) 
  
  # reshape data matrices for easier plotting with ggplot. 
  # NB: the names of all variables in the reshaped data frames should match, 
  # otherwise there will be a ggplot scale error (the average profiles data frame
  # has less rows than the full profiles data frame)! 
  div.prof.long <- melt(div.profiles)
  names(div.prof.long) <- c("q", "replicate", "value")
  
  aver.prof.long <- melt(aver.profiles)
  names(aver.prof.long) <- c("q", "replicate", "value")
  
  # add in the stations to serve as grouping factor
  
  # the number of values for each station corresponds to the number of values of
  # q (= nrows of data matrix) * number of replicates per station (assumes that
  # there is the same number of replicates per station).
  div.prof.long$station <- rep(stations, 
                               each = nrow(div.profiles) * ncol(div.profiles)/length(stations))
  
  aver.prof.long$station <- rep(stations, 
                                each = nrow(aver.profiles) * ncol(aver.profiles)/length(stations))
  
  # reorder the new factor's levels according to the order of the input factor 
  div.prof.long$station <- reorder.factor(div.prof.long$station,
                                          new.order = stations)
  aver.prof.long$station <- reorder.factor(aver.prof.long$station,
                                           new.order = stations)
  
  # plot 
  p <- ggplot(data = div.prof.long, aes(x = q,
                                        y = value, 
                                        group = replicate)) + 
    
    # plot the diversity profiles for each replicate 
    geom_line(colour = col.profiles[1], show.legend = FALSE) +
    
    # plot the average diversity profile for each station
    geom_line(data = aver.prof.long, lwd = 0.8, colour = col.profiles[2]) +
    facet_wrap(~station) +
    labs(x = "Sensitivity parameter q", y = "Diversity") +
    theme_bw()
  
  return(p)
}