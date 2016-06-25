plot_contrib_tax_groups <- function(tax.gr.props, by.years = FALSE) {
  ## prepares and plots contributon of taxonomic groups to the abundance/
  ## biomass (as proportions). 
  
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
        labs(x = "", y = "Proportion") + 
        theme_bw()
  
  # facet by year, if specified at input
  if(by.years){
    p <- p + facet_wrap(~years, nrow = 2)
  }
  
  return(p)
  
}

p <- ggplot(data = y, aes(x = stations, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position="stack") +
  #custom.col.scale +
  scale_fill_viridis(discrete = TRUE, name = "") +
  labs(x = "", y = "Proportion") + 
  theme_bw()
