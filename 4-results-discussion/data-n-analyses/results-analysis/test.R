plot_contrib_tax_groups <- function(tax.gr.props, by.years = FALSE) {
  ## prepares and plots contributon of taxonomic groups to the abundance/
  ## biomass (as proportions). 
  
  library(ggplot2)
  library(plyr)
  library(reshape2)
  
  
  if(by.years) {
    # summarize data by stations and years, if specified at input 
    summary.tax.props <- ddply(tax.gr.props, .(stations, years), colwise(mean, .cols = is.numeric))  
    
  } else {
    # summarize data by stations only
    summary.tax.props <- ddply(tax.gr.props, .(stations), colwise(mean, .cols = is.numeric))  
  }
  
  # reshape the data frame in long format
  tax.props.melted <- melt(summary.tax.props)  
  
  # define custom color scheme based on the colours specified at input 
  custom.cols <- col.curves
  names(custom.cols) <- levels(curves.melted$variable)
  custom.col.scale <- scale_color_manual(name = "", values = custom.cols)
  
  # plot
  p <- ggplot(data = tax.props.melted, aes(x = stations, y = value, fill = variable)) + 
        geom_bar(stat="identity", position="stack") +
        custom.col.scale +
        labs(x = "Stations", y = "Proportion", colour = "") + 
        theme_bw()
  
  # facet by year, if specified at input
  if(by.years){
    p <- p + facet_wrap(~years)  
  }
  
  
  return(p)
  
}