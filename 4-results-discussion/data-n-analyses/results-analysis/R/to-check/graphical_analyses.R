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
abc <- function (abnd.data, biomass.data, abnd.val, biomass.val, ...) {
  ## calculates the cumulative abundance and biomass percentages for plotting ABC 
  ## curves. 
  ## Arguments: abnd.data, biomass.data - long tibbles of zoobenthic abundance and biomass.
  ##            abnd.val, biomass.val - name of the column containing the numeric 
  ## values for the abundance/biomass
  ##            ... - any grouping variables. NB should be THE SAME in both abundance and biomass dfs
  ## Dependencies: tidyverse packages, vegan
  
  
  # get the columns containing the actual values of the abundance & the biomass
  abnd.col <- enquo(abnd.val)
  biomass.col <- enquo(biomass.val)
  
  # get the grouping variables
  group.vars = quos(...)
  
  # sort the abundance and biomass data in decreasing order, KEEPING THE OTHER 
  # COLUMNS IN ORDER BY STATION, YEAR, MONTH & REPLICATE - OTHERWISE THE ABUNDANCE 
  # & BIOMASS WILL NEVER MATCH EXACTLY!
  abnd <- abnd.data %>% 
    group_by(!!!group.vars) %>%
    arrange(desc(!!abnd.col), .by_group = TRUE) 
  
  biomass <- biomass.data %>% 
    group_by(!!!group.vars) %>%
    arrange(desc(!!biomass.col), .by_group = TRUE) 

  # check if the abundance and biomass are imported correctly - they should 
  # contain the same number of elements unless there is an error in the input 
  # data - e.g. a species for which there is only an abundance or only a 
  # biomass.
  if (nrow(abnd) != nrow(biomass)) {
    stop("The abundance and biomass datasets contain different numbers of species! Check your input data!")
  }
  
  # calculate the number of species S in each replicate (will be needed later 
  # for the W statistic). Shouldn't matter which dataset is used - abundance or biomass
  s <- abnd %>%
    summarize(s = vegan::specnumber(!!abnd.col))

  # calculate the cumulative percent abundance for each species in each replicate
  abnd <- abnd %>%
    mutate(perc_abnd = cumsum(!!abnd.col) / sum(!!abnd.col) * 100)
  
  biomass <- biomass %>%
    mutate(perc_biomass = cumsum(!!biomass.col) / sum(!!biomass.col) * 100)

  # make new tibble with the results - by binding, otherwise will get mixed up because of multiple 
  # identical combinations. MAKE SURE THE OBSERVATIONS ARE IN THE RIGHT ORDER!

  if(!identical(abnd %>% select(!!!group.vars), 
                biomass %>% select(!!!group.vars))) {
    
    stop("Observations in abundance & biomass datasets not in the same order!")
  
  } else {
    
    abnd.sub <- abnd %>% select(!!!group.vars, perc_abnd)
    biomass.sub <- biomass %>% select(perc_biomass)
    
    abc <- bind_cols(abnd.sub, biomass.sub)
    
    abc <- abc %>% 
      select(!!!group.vars, perc_abnd, perc_biomass) # the grouping columns from the biomass tibble are added by default, but are not really needed
    
  }

  ## calculate the W statistic (= quantification of the distance between the abundance
  ## & the biomass curves)
  w <- abc %>%
    summarize(sum_ba = sum((perc_biomass - perc_abnd)))
  
  w <- left_join(w, s) %>%
    mutate(w = sum_ba / (50 * (s - 1))) %>%
    select(!!!group.vars, w)
  

  return(list(abc, w))
}



#### STILL HAVE TO FIX: #####

partial_dominance_curves <- function(abnd.data, biomass.data) {
  ## calculates partial dominance curves from a macrobenthic abundance-biomass
  ## dataset. Arguments: abnd.data, biomass.data - vectors/data frames of abundance
  ## and biomass, respectively.
  
  # get the abundance and biomass data and sort them in decreasing order 
  abundance <- sort(abnd.data, decreasing = TRUE)
  biomass <- sort(biomass.data, decreasing = TRUE)
  
  # check if the abundance and biomass vectors contain the same number of elements
  # (i.e., if there is an error in the input data - a species for which only an
  # abundance or only a biomass value is given)
  if (sum(abundance != 0) != sum(biomass != 0)) {
    stop("The abundance and biomass vectors have different lengths! Check your input data!")
  }
  
  # get rid of double-0 - species absent from the sample, which would lead to a division by 0
  abundance <- abundance[abundance > 0]
  biomass <- biomass[biomass > 0]
  
  # make a vector to hold partial percent abundances 
  perc.abnd <- NA
  # calculate the partial abundances, successively removing each species and calculating 
  # over the rest 
  for (i in 1:length(abundance)) {
    
    if (sum(abundance[i:length(abundance)]) == 0){
      perc.abnd[i] <- NA
    }
    else {
      perc.abnd[i] <- (abundance[i] / sum(abundance[i:length(abundance)])) * 100  
    }    
  }
  
  # same procedure for the biomass
  perc.bio <- NA
  for (i in 1:length(biomass)) {
    if (sum(biomass[i:length(biomass)]) == 0){
      perc.bio[i] <- NA
    }
    else {
      perc.bio[i] <- (biomass[i] / sum(biomass[i:length(biomass)])) * 100  
    }      
  }
  
  # bind the partial abundance and biomass into a new data frame, without matching
  # the species names to its specific values for partial abundance & biomass.
  abc.partial <- as.data.frame(cbind(perc.abnd, perc.bio))
  colnames(abc.partial) <- c("partial.abnd", "partial.biomass")
  
  return(abc.partial)
}
