###### Example graphs for M&M - k-dominance and ABC curves.

## define working subdirectories
data.dir <- "data"
save.dir <- "output"
functions.dir <- "R"
figures.dir <- "figs"

###### k-dominance curves. Data: Zettler, 2005: Baltic Sea macrozoobenthos (OBIS).
###### http://www.iobis.org/mapper/?dataset=254

## import & clean, keeping only some of the data (only want it for illustrative purposes).
baltic.sea.zoo <- read.csv(file.path(data.dir, "dominance_curves", "Zettler et al_2005_OBIS_Baltic_macrozoobenthos_data.csv"), 
                           header = TRUE, stringsAsFactors = FALSE)

str(baltic.sea.zoo)

baltic.sub <- baltic.sea.zoo[, c("yearcollected", "monthcollected", "tname", 
                                 "locality", "observedindividualcount")]
str(baltic.sub)

baltic.sub[baltic.sub$yearcollected == 2005,] # only 2 stations sampled - not enough..
unique(baltic.sub[baltic.sub$yearcollected == 2002,]$locality) # 7 stations - prob. too many..
unique(baltic.sub[baltic.sub$yearcollected == 1990,]$locality) # 5 stations - maybe ok

## try plotting k-dominane curves first for 1990, then 2002, to decide which looks better.
baltic.1990 <- baltic.sub[baltic.sub$yearcollected == 1990,]
baltic.2002 <- baltic.sub[baltic.sub$yearcollected == 2002,]

# subset again, to leave only species, stations and counts
baltic.1990.sub <- baltic.1990[, c("tname", "locality", "observedindividualcount")] 
baltic.2002.sub <- baltic.2002[, c("tname", "locality", "observedindividualcount")]

# make localities into factors
baltic.1990.sub$locality <- as.factor(baltic.1990.sub$locality)
baltic.2002.sub$locality <- as.factor(baltic.2002.sub$locality)

# assign ranks to species 
library(plyr)
baltic.1990.sub.prop <- ddply(baltic.1990.sub, .(locality), function(x) {
  # calculate relative abundance of each species (proportion of the total abundance)
  x$prop = x$observedindividualcount / sum(x$observedindividualcount)
  # assign rank to species
  x <- x[order(x$prop, decreasing = TRUE),]
  x$rank <- 1:nrow(x)
  
  # calculate the cumulative proportion (easier plotting)
  x$cumprop[1] <- x$prop[1]
  for (i in 2:nrow(x)) {
    x$cumprop[i] <- x$cumprop[i-1] + x$prop[i]
  }
   
  return(x)
})


baltic.2002.sub.prop <- ddply(baltic.2002.sub, .(locality), function(x) {
  # calculate relative abundance of each species (proportion of the total abundance)
  x$prop = x$observedindividualcount / sum(x$observedindividualcount)
  # assign rank to species
  x <- x[order(x$prop, decreasing = TRUE),]
  x$rank <- 1:nrow(x)
  
  # calculate the cumulative proportion (easier plotting)
  x$cumprop[1] <- x$prop[1]
  for (i in 2:nrow(x)) {
    x$cumprop[i] <- x$cumprop[i-1] + x$prop[i]
  }
  
  return(x)
})


# plot dominance curves
library(ggplot2)
library(viridis)

baltic.1990.plot <- ggplot(baltic.1990.sub.prop, 
                           mapping = aes(x = rank, y = cumprop * 100,
                                         colour = locality, shape = locality)) + 
                      geom_line() +
                      geom_point() +
                      scale_x_log10() +
                      scale_color_viridis(discrete = TRUE, option = "viridis") + 
                      ylab("Cumulative abundance (%)") +
                      xlab("Species rank (log)") +
                      theme_bw() + 
                      theme(legend.position = "none")

# yells about too many series for the shapes palette (max 6, here - 7), but works
baltic.2002.plot <- ggplot(baltic.2002.sub.prop, 
                           mapping = aes(x = rank, y = cumprop * 100,
                                         colour = locality, shape = locality)) + 
  geom_line() +
  geom_point() + 
  scale_x_log10() +
  scale_color_viridis(discrete = TRUE, option = "viridis") + 
  ylab("Cumulative abundance (%)") +
  xlab("Species rank (log)") +
  theme_bw() + 
  theme(legend.position = "none")


ggsave(file.path(figures.dir, "kdom-baltic-1990.png"), baltic.1990.plot, 
       dpi = 500, width = 13, units = "cm")
ggsave(file.path(figures.dir, "kdom-baltic-2002.png"), baltic.2002.plot, 
       dpi = 500, width = 13, units = "cm")


## cleanup
rm(baltic.1990, baltic.1990.sub, baltic.1990.sub.prop, baltic.2002, baltic.2002.sub,
   baltic.2002.sub.prop, baltic.sea.zoo, baltic.sub)

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
