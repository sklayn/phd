plot_div_profiles <- function(div.profiles, stations, one.panel = FALSE) {
  ## plots the diversity profiles by sampling station.
  ## Arguments: div.profiles - matrix of diversity profiles (in wide format). 
  ##              Each profile is a column (rows = q values, columns = 
  ##              samples/stations).
  ##            stations - vector of station names. The order of its levels is 
  ##              assumed to be the correct one, to be used in the plot.
  ##            one.panel - logical, should the profiles for each station be 
  ##              placed in a single panel (for example when average profiles 
  ##              only), or in separate panels by station? Default - separate 
  ##              panel for each station.
  ## Output: ggplot object.
  
  # import the necessary libraries
  library(ggplot2)
  library(plyr)
  library(reshape2)  

  # reshape data matrix for easier plotting with ggplot
  div.prof.long <- melt(div.profiles)
  names(div.prof.long) <- c("q", "replicate", "value")
  
  # add in the stations to serve as grouping factor

  nq <- nrow(div.profiles)                      # number of q values
  nrep <- ncol(div.profiles)/length(stations)   # number of replicates/station
  
  div.prof.long$station <- rep(stations.sand, each = nq * nrep)
  
  # reorder the new factor's levels according to the order of the input factor 
  div.prof.long$station <- reorder.factor(div.prof.long$station,
                                          new.order = stations.sand)
  
  if (one.panel) {
    # plot the profiles in one panel (e.g., when one average profile/station)  
    ggplot(data = div.prof.long, aes(x = q, 
                                     y = value, 
                                     group = replicate, 
                                     colour = station)) + 
      geom_line() +
      labs(x = "Sensitivity parameter q", y = "Diversity", colour = "Stations") +
      theme_bw()  
    
  } else {
    # plot the profiles in separate panels by station  
    ggplot(data = div.prof.long, aes(x = q, 
                                     y = value, 
                                     group = replicate, 
                                     colour = station)) + 
      geom_line(show.legend = FALSE) +
      facet_wrap(~station) +
      labs(x = "Sensitivity parameter q", y = "Diversity") +
      theme_bw()  
  }

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
  div.prof.long$station <- rep(stations.sand, 
                               each = nrow(div.profiles) * ncol(div.profiles)/length(stations))
  
  aver.prof.long$station <- rep(stations.sand, 
                                each = nrow(aver.profiles) * ncol(aver.profiles)/length(stations))

  # reorder the new factor's levels according to the order of the input factor 
  div.prof.long$station <- reorder.factor(div.prof.long$station,
                                          new.order = stations.sand)
  aver.prof.long$station <- reorder.factor(aver.prof.long$station,
                                           new.order = stations.sand)
  
  # plot 
  ggplot(data = div.prof.long, aes(x = q,
                                   y = value, 
                                   group = replicate)) + 
    
    # plot the diversity profiles for each replicate 
    geom_line(colour = col.profiles[1], show.legend = FALSE) +
    
    # plot the average diversity profile for each station
    geom_line(data = aver.prof.long, lwd = 0.8, colour = col.profiles[2]) +
    facet_wrap(~station) +
    labs(x = "Sensitivity parameter q", y = "Diversity") +
    theme_bw()

}