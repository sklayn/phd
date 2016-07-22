summary_zoo_params <- function(num.data, factors, param = c("abnd", "biomass", "nb.sp")) {
  ## helper function to summarize (quick and dirty) the main parameters of a 
  ## zoobenthic dataset: average abundance, biomass or number of species, 
  ## min/max of those, and which samples correspond to the min/max values.
  ## 
  ## Arguments: num.data - samples x species data frame of abundance/biomass
  ##            factors - factors data frame, where the rows are in the same 
  ##             order as the data frame of numeric data - containing station id, 
  ##             years/dates, etc
  ##            param - which parameter should be summarized. 
  
  library(vegan)
  
  if(param == "abnd" | param == "biomass") {
    # calculate mean and standard deviation of the parameter
    mean.val <- mean(rowSums(num.data))
    sd.val <- sd(rowSums(num.data))
    
    # find the minimum value, and look it up in the factors data frame 
    min.val <- min(rowSums(num.data))
    min.id <- factors[which.min(rowSums(num.data)), ]
    
    # find the maximum value, and look it up
    max.val <- max(rowSums(num.data))
    max.id <- factors[which.max(rowSums(num.data)), ]
    
  } else {
    # use specnumber from package vegan if summarizing the number of species; 
    # same procedure as for the abundance/biomass
    mean.val <- mean(specnumber(num.data))
    sd.val <- sd(specnumber(num.data))
    
    min.val <- min(specnumber(num.data))
    min.id <- factors[which.min(specnumber(num.data)), ]
    
    max.val <- max(specnumber(num.data))
    max.id <- factors[which.max(specnumber(num.data)), ]
  }
  
  # put the mean and sd values in a data frame (cleaner)
  mean.sd <- data.frame(mean = mean.val, st.dev = sd.val)
  
  # put the minimum-maximum values in a list of their own
  min.all <- list("min.value" = min.val, "min.id" = min.id)
  max.all <- list("max.value" = max.val, "max.id" = max.id)
  
  # combine and return all of these as a list  
  summary.list <- list("mean.sd" = mean.sd, 
                       "min" = min.all,
                       "max" = max.all) 
  
  return(summary.list)
}