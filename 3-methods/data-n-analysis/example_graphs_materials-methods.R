###### Example graphs for M&M - k-dominance and ABC curves.

## define working subdirectories
data.dir <- "data"
save.dir <- "output"
functions.dir <- "R"
figures.dir <- "figs"

###### k-dominance curves. Data: 





###### ABC curves. Data: Loch Linhe abundance & biomass (example data from PRIMER 6 
###### manual). Since this is only for illustrative purposes, 3 partiuclar years 
###### are chosen: 1965 (normal), 1970 (moderately polluted with organic matter 
###### from pulp mill effluent), and 1972 (very polluted). 
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
# source(file.path("/media/eti/My Passport/rabota/documents/phd/thesis/4-results-discussion/data-n-analyses/results-analysis",
#                  functions.dir, "plot_dom_curves.R"))

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
  
# get rid of the third column, BiAi - not needed here 
linhe.abc.for.plot <- lapply(linhe.abc.for.plot, "[", -3)

## reshape and plot the curves - have to fix/get rid of the function later; 
## it's horrible
library(ggplot2)
library(reshape2)
library(plyr)

linhe.abc.for.plot <- lapply(linhe.abc.for.plot, function(x) {
  names(x) <- c("abundance", "biomass")
  x$rank <- 1:nrow(x)
  return(x)
})

## convert to data frame & reshape 
linhe.abc.df <- ldply(linhe.abc.for.plot)
linhe.abc.df.long <- melt(linhe.abc.df, id.vars = c(".id", "rank"))

## plot
linhe.abc.plot <- ggplot(linhe.abc.df.long) + 
    geom_line(aes(x = rank, y = value, colour = variable)) +
    scale_color_manual(values = c("skyblue", "orange")) +
    scale_x_log10() +
    labs(x = "Species rank", y = "Cumulative %", colour = "") +
    facet_wrap(~.id, ncol = 3) +
    theme_bw()

ggsave(file.path(figures.dir, "loch-linhe-abc-curves.png"), 
       plot = linhe.abc.plot, 
       width = 16, height = 8, units = "cm", dpi = 500)  

## cleanup
rm(linhe.abc, linhe.abc.df, linhe.abc.df.long, linhe.abc.for.plot, linhe.1965.abc, 
   linhe.abnd, linhe.abnd.fin, linhe.abnd.sub, linhe.biomass, linhe.biomass.fin, 
   linhe.biomass.sub)
