plot_dom_curves <- function(dom.curves, trasf.y = FALSE, col.curves = c("skyblue", "orange")) {
  ## plots ABC or partial dominance curves for stations/replicates. 
  ## Arguments: dom.curves - list of curves to be plotted (each element = 
  ##              station/replicate).
  ##            trasf - should the y axis should be transformed using a (modified)
  ##              logistic transformation for a better visualization of the curves?
  ##            col.curves - custom colours for the abundance and biomass curves,
  ##              respectively.
  
  library(ggplot2)
  library(reshape)
  library(scales)

  # the input is still a list, so convert to data frame and drop the third column
  # if there is one (as in the case of ABC)
  if(class(dom.curves) == "list") {
    curves <- as.data.frame(dom.curves)[-3]  
  }
  
  names(curves) <- c("abundance", "biomass")
  
  # add a variable for the species rank
  curves$rank <- 1:nrow(curves)
  
  # reshape the data frame in long format
  curves.melted <- melt(curves, id.vars = "rank")  
  
  # make a custom colour scale based on the colours specified at input 
  custom.cols <- col.curves
  names(custom.cols) <- levels(curves.melted$variable)
  custom.col.scale <- scale_color_manual(name = "", values = custom.cols)
  
  if(trasf.y) {
    # plot using the modified y scale
    ggplot(data = curves.melted, aes(x = rank, y = value, colour = variable)) + 
      geom_line() +
      custom.col.scale +
      scale_x_log10() +
      coord_trans(y = "modif_logistic") +
      labs(x = "Species rank", y = "Cumulative %", colour = "") +
      theme_bw()
    
  } else {
    # plot without transforming the y axis
    ggplot(data = curves.melted, aes(x = rank, y = value, colour = variable)) + 
      geom_line() +
      custom.col.scale +
      scale_x_log10() +
      labs(x = "Species rank", y = "Cumulative %", colour = "") +
      theme_bw()  
  }
  
}


plot_dom_curves_facets <- function(dom.curves, stations, trasf.y = FALSE, col.curves = c("skyblue", "orange")) {
  ## plots ABC or partial dominance curves for stations/replicates. 
  ## Arguments: dom.curves - list of k-dominance curves (abundance/biomass) for 
  ##             stations/replicates.
  ##            stations - vector of station names for the current dataset; length
  ##              should match the number of elements in dom.curves (= number of 
  ##              stations/replicates), and they shouold be in the same order, too.
  ##            trasf - should the x axis should be transformed using a modified 
  ##              logistic transformation (for better visualization of the curves)?
  ##            col.curves - custom colours for the abundance and biomass curves,
  ##              respectively.
  
  library(ggplot2)
  library(reshape)
  library(scales)
  library(plyr)
  
  # reshape the input (list of lists) into one big data frame
  curves <- ldply(dom.curves, data.frame)
  
  if(ncol(curves) > 3) {
    # if there are more than 3 columns after reshaping the list into data frame 
    # (happens with ABC, which gets columns for Bi-Ai and W, too) - drop the 
    # unnecessary columns. The first 3 are all we need (= id, abundance and 
    # biomass).
    curves <- curves[1:3]  
  }
  
  names(curves) <- c("id", "abundance", "biomass")
  curves$id <- as.factor(curves$id)
  
  # add a variable for species rank - will serve as x axis in the plots
  curves <- ddply(curves, "id", transform, rank = 1:length(id))
  
  # reshape the data frame in long format
  curves.melted <- melt(curves, id.vars = c("id", "rank"))  
  
  # rename the id variable to correspond to the station names supplied at input
  levels(curves.melted$id) <- stations
  
  # make a custom colour scale based on the colours specified at input 
  custom.cols <- col.curves
  names(custom.cols) <- levels(curves.melted$variable)
  custom.col.scale <- scale_color_manual(name = "", values = custom.cols)
  
  if(trasf.y) {
    # plot using the modified y scale
    ggplot(data = curves.melted, aes(x = rank, y = value, colour = variable)) + 
      geom_line() +
      custom.col.scale +
      scale_x_log10() +
      coord_trans(y = "modif_logistic") +
      labs(x = "Species rank (log)", y = "Cumulative %", colour = "") +
      facet_wrap(~id)+
      theme_bw()
    
  } else {
    # plot without transforming the y axis
    ggplot(data = curves.melted, aes(x = rank, y = value, colour = variable)) + 
      geom_line() +
      custom.col.scale +
      scale_x_log10() +
      labs(x = "Species rank (log)", y = "Cumulative %", colour = "") +
      facet_wrap(~id) +
      theme_bw()  
  }
  
}
