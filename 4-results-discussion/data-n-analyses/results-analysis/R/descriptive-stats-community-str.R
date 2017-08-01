##### Descriptive community statistics #####

summary_zoo_params <- function(num.data, metadata, parameter = c("abnd", "biomass", "nb.sp")) {
  ## Summarizes (quick and dirty) the main parameters of a zoobenthic dataset: 
  ##  average abundance, biomass or number of species, min/max of those, and 
  ##  which samples correspond to the min/max values.
  ## 
  ## Arguments: num.data - samples x species data frame of abundance/biomass
  ##            metadata - data frame of descriptive metadata (station id, 
  ##             month, year, etc.), where the rows correspond to those 
  ##             in the abundance/biomass data frame.
  ##            param - which parameter(s) should be summarized.
  ## Output: list of the summary data.
  ## Dependencies: vegan (for number of species).
  
  library(vegan)
  
  if(parameter == "abnd" | parameter == "biomass") {
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
  summary.list <- list("mean_sd" = mean.sd, 
                       "min" = min.all,
                       "max" = max.all) 
  
  return(summary.list)
}



##### Contribution of taxonomic groups to the abundance/biomass #####
tax_group_contribution <- function(abnd.data, tax.groups) {
  ## Calculates the proportional contribution of each of a set of larger taxonomic
  ##  groups to the overall abundance/biomass in the stations' communities.
  ## Arguments: abnd.data - data frame of numeric abundances/biomasses for all
  ##              taxa present in the current dataset (stations x species)
  ##            tax.data - vector of taxonomic group labels of the same length 
  ##              as the number of non-0 (i.e., present) species in the current
  ##              dataset
  ##           
  ## Output: data frame of proportions of each taxonomic group in each station. 
  ##  For each station, the sum of all taxonomic groups' proportions (= rows) 
  ##  should be 1.
  ## Dependencies: plyr 
  
  library(plyr)
  
  # calculate each species' proportion in each sample 
  prop.abnd <- abnd.data / rowSums(abnd.data)
  
  # transpose the proportions data frame and add the factor column containing each
  # species' taxonomic group
  prop.abnd <- as.data.frame(t(prop.abnd))
  names(prop.abnd) <- 1:ncol(prop.abnd)
  
  tax.abnd <- cbind(tax.group = tax.groups, prop.abnd)
  
  # aggregate data by taxonomic group (= sum of proportions of all sp in that group)
  summary.tax.abnd <- ddply(tax.abnd, .(tax.group), colwise(sum, .cols = is.numeric))
  
  # transpose to a final data frame of stations x taxonomic groups 
  summary.tax.abnd <- setNames(as.data.frame(t(summary.tax.abnd[, -1])), summary.tax.abnd[, 1])
  
  return(summary.tax.abnd)
}


plot_tax_group_contribution <- function(tax.gr.props, by.years = FALSE) {
  ## Prepares and plots in ggplot2 the contributon of taxonomic groups to the 
  ## abundance/biomass (exressed as proportions) of the communities. 
  ## Arguments: tax.gr.props - data frame (in wide format) of taxonomic groups'
  ##              proportions; stations x taxonomic groups, with additional 
  ##              identifying factors (expected: stations and years; all others
  ##              are ignored by default).
  ##            by.years - should the final plot be faceted by year?
  ## Returns ggplot2 object.
  
  library(ggplot2)
  library(plyr)
  library(reshape2)
  library(viridis)
  
  if(by.years) {
    # summarize data by stations and years, if specified at input 
    summary.tax.props <- ddply(tax.gr.props, .(stations, years), colwise(mean, .cols = is.numeric))  
    
  } else {
    # summarize data by stations only
    summary.tax.props <- ddply(tax.gr.props, .(stations), colwise(mean, .cols = is.numeric))  
  }
  
  # reshape the data frame in long format
  tax.props.melted <- melt(summary.tax.props)  
  
  # plot proportions as stacked bars
  p <- ggplot(data = tax.props.melted, aes(x = stations, y = value, fill = variable)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_viridis(discrete = TRUE, name = "") +
    labs(x = "Station", y = "Proportion") + 
    theme_bw() + 
    # modify the size of x axis and legend labels 
    theme(axis.text.x = element_text(size = rel(1.3)),
          axis.text.y = element_text(size = rel(1.3)),
          axis.title.x = element_text(size = rel(1.2)),
          legend.text = element_text(size = rel(1.1)))
  
  # facet by year, if specified at input
  if(by.years){
    p <- p + facet_wrap(~years, nrow = 2) + 
      theme(strip.text.x = element_text(size = rel(1.3)))
  }
  
  return(p)
}
