###### Example graphs for M&M - k-dominance and ABC curves.

## define working subdirectories
data.dir <- "data"
save.dir <- "output"
functions.dir <- "R"
figures.dir <- "figs"

###### Data: Loch Linhe abundance & biomass (example data from PRIMER 6 manual).
###### Since this is only for illustrative purposes, 3 partiuclar years are 
###### chosen: 1965 (normal), 1970 (moderately polluted with organic matter from
###### pulp mill effluent), and 1972 (very polluted). 
linhe.abnd <- read.delim(file.path(data.dir, "dominance_curves", 
                                   "lnma.txt"), 
                         sep = "\t", header = TRUE)

linhe.biomass <- read.delim(file.path(data.dir, "dominance_curves", 
                                      "lnmb.txt"), 
                            sep = "\t", header = TRUE)

## fix col.names
names(linhe.abnd) <- gsub(pattern = "X", replacement = "", names(linhe.abnd))
names(linhe.biomass) <- gsub(pattern = "X", replacement = "", names(linhe.biomass))

## subset to only 1965, 1970 & 1972
linhe.abnd.sub <- linhe.abnd[, c("1965", "1970", "1972")] 
linhe.biomass.sub <- linhe.biomass[, c("1965", "1970", "1972")]


#### source ABC calculation/plotting functions (from results/discussion dir)
source(file.path("/media/eti/My Passport/rabota/documents/phd/thesis/4-results-discussion/data-n-analyses/results-analysis", 
                 functions.dir, "k-dom_curves.R"))
source(file.path("/media/eti/My Passport/rabota/documents/phd/thesis/4-results-discussion/data-n-analyses/results-analysis", 
                 functions.dir, "plot_dom_curves.R"))

#### THIS VARIANT OF THE ABC FUNCTION MIGHT NEED TO BE FIXED - DEMANDS BOTH ABND AND 
#### BIOMASS FOR EACH SPECIES IN A SAMPLE, AND THE LINHE DATA DOESN'T ALWATYS HAVE THAT.
#### SO, WHERE THE BIOMASS WAS 0, IT WAS CHANGED (MANUALLY) TO 0.0001, FOR SPEED
linhe.abnd.fin <- read.csv(file.path(data.dir, "dominance_curves", 
                                   "linhe.abnd.sub.csv"), 
                          header = TRUE)

linhe.biomass.fin <- read.csv(file.path(data.dir, "dominance_curves", 
                                        "linhe.bio.sub.csv"), 
                              header = TRUE)

## fix col.names again..
names(linhe.abnd.fin) <- gsub(pattern = "X", replacement = "", names(linhe.abnd.fin))
names(linhe.biomass.fin) <- gsub(pattern = "X", replacement = "", names(linhe.biomass.fin))


linhe.1965.abc <- abc(linhe.abnd.fin$`1965`, linhe.biomass.fin$`1965`)

linhe.abc <- mapply(abc, linhe.abnd.fin, linhe.biomass.fin, SIMPLIFY = FALSE)
linhe.abc.for.plot <- lapply(linhe.abc, function(x) x[[1]]) 
  
  
## quickly fix the stupid plotting function here, because no more time & nerves
plot_dom_curves2 <- function(dom.curves, trasf.y = FALSE, col.curves = c("skyblue", "orange")) {
  ## plots ABC or partial dominance curves for stations/replicates. 
  ## Arguments: dom.curves - list of curves to be plotted (each element = 
  ##              station/replicate).
  ##            trasf - should the y axis should be transformed using a (modified)
  ##              logistic transformation for a better visualization of the curves?
  ##            col.curves - custom colours for the abundance and biomass curves,
  ##              respectively.
  
  library(ggplot2)
  library(reshape2)
  library(scales)
  library(plyr)
  
  # the input is still a list, so convert to data frame and drop the third column
  # if there is one (as in the case of ABC)
  curves <- dom.curves[-3]  
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
  
abc.plots <- lapply(linhe.abc.for.plot, plot_dom_curves2)  

library(gridExtra)
grid.arrange(abc.plots[[1]], abc.plots[[2]], abc.plots[[3]])