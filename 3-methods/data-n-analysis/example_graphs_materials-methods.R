###### Example graphs for M&M - statistical analyses ######

## define working subdirectories
data.dir <- "data"
save.dir <- "output"
functions.dir <- "R"
figures.dir <- "figs"

## import libraries
library(ggplot2)
library(mvabund)
library(plyr)
library(reshape2) # deprecated? see tidyr
library(tidyr)
library(viridis)

###### k-dominance curves ##### 
#### Data: Zettler, 2005: Baltic Sea macrozoobenthos (OBIS).
#### http://www.iobis.org/mapper/?dataset=254

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

###### ABC curves ######  
#### Data: Loch Linhe abundance & biomass (example data from PRIMER 6 manual). 
#### Since this is only for illustrative purposes, 3 partiuclar years are chosen: 
#### 1965 (normal), 1970 (moderately polluted with organic matter from pulp mill 
#### effluent), and 1972 (very polluted). 
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


##### Mean-variance relationships in multivariate abundance datasets #####
#### Data: Zettler, 2005: Baltic Sea macrozoobenthos (OBIS).
#### http://www.iobis.org/mapper/?dataset=254

## import & clean, keeping only some of the data (only want it for illustrative purposes).
baltic.sea.zoo <- read.csv(file.path(data.dir, "dominance_curves", "Zettler et al_2005_OBIS_Baltic_macrozoobenthos_data.csv"), 
                           header = TRUE, stringsAsFactors = FALSE)

str(baltic.sea.zoo)

baltic.sub <- baltic.sea.zoo[, c("yearcollected", "monthcollected", "tname", 
                                 "locality", "observedindividualcount")]
str(baltic.sub)

unique(baltic.sub[baltic.sub$yearcollected == 2002,]$locality) # 7 stations sampled - will do for illustraton.
unique(baltic.sub[baltic.sub$yearcollected == 2003,]$locality) # 9 stations sampled - will do for illustraton.

## get only 2003, and only the taxon, station and count
baltic.sub.2003 <- baltic.sub[baltic.sub$yearcollected == 2003, c("tname", "locality", "observedindividualcount")]

## fill in missing combinations (add 0-count taxa to whatever station they're missing from) 
## (package tidyr)
baltic.sub.2003 <- tidyr::complete(baltic.sub.2003,
                            tname, locality,
                            fill = list(observedindividualcount = 0))

## remove Moni-018; it's causing trouble with duplicate values that I don't have time to diagnose now
baltic.sub.2003.fixed <- baltic.sub.2003[-which(baltic.sub.2003$locality == "Moni-018"), ]

# calculate mean and variance for each taxon, then plot in ggplot
baltic.sub.2003.fixed.gg <- plyr::ddply(baltic.sub.2003.fixed, 
                                  .variables = "tname", 
                                  summarize, 
                                  mean = mean(observedindividualcount), 
                                  variance = var(observedindividualcount))

meanvar.plot.baltic.2003 <- ggplot(baltic.sub.2003.fixed.gg) + 
  geom_point(aes(x = mean + 1, y = variance + 1)) + # 1 added to avoid log(0) = -Inf
  scale_x_continuous(trans = "log10", name = "Mean (log)") +
  scale_y_continuous(trans = "log10", name = "Variance (log)") +
  theme_bw() + 
  ggtitle("Untransformed counts") + 
  theme(plot.title = element_text(hjust = 0.5))

## try log(x + 1) transforming the taxon counts before plotting
baltic.sub.2003.fixed.gg.log <- plyr::ddply(baltic.sub.2003.fixed, 
                                  .variables = "tname", 
                                  summarize, 
                                  mean = mean(log(observedindividualcount + 1)), 
                                  variance = var(log(observedindividualcount + 1)))

meanvar.plot.baltic.2003.log <- ggplot(baltic.sub.2003.fixed.gg.log) + 
  geom_point(aes(x = mean + 1, y = variance + 1)) + # 1 added to avoid log(0) = -Inf
  scale_x_log10(name = "Mean (log)", limits = c(1, 1000)) + # same scale as untransformed counts plot
  scale_y_log10(name = "Variance (log)", limits = c(1, 10000000)) + # idem x axis
  theme_bw() +
  ggtitle("log(x + 1) counts") +
  theme(plot.title = element_text(hjust = 0.5))
    


meanvar.plots <- gridExtra::grid.arrange(meanvar.plot.baltic.2003,
                              meanvar.plot.baltic.2003.log,
                              ncol = 2)

gridExtra::grid.arrange(meanvar.plot.baltic.2003, 
             meanvar.plot.baltic.2003.log, 
             ncol = 2)

ggsave(file.path(figures.dir, "mean-variance-ex-plot.png"), 
       meanvar.plots,
       width = 16.5, units = "cm", 
       dpi = 500)


### plot the mean-variance relationship (log scale) - using package mvabund
# first, extract only the abundances, transpose and make into mvabund matrix
library(reshape2) 
baltic.sub.2003.fin <- reshape2::dcast(baltic.sub.2003.fixed, 
                             formula = locality ~ tname, 
                             value.var = "observedindividualcount")
# plot
library(mvabund)
meanvar.plot(mvabund(baltic.sub.2003.fin[, -1]),
             xlab = "Mean (log scale)", 
             ylab = "Variance (log scale)")

# no log-transform of axes
meanvar.plot(mvabund(baltic.sub.2003.fin[, -1]), log = "",
             xlab = "Mean", ylab = "Variance")

## transform data (log(x + 1)) & plot again
baltic.sub.2003.fin.log <- apply(baltic.sub.2003.fin[, -1], 
                                 2, 
                                 function(x) log(x + 1))

meanvar.plot(baltic.sub.2003.fin.log, 
             xlab = "Mean (log scale)", 
             ylab = "Variance (log scale)")

### just use the data/examples that come with mvabund
## Load the tikus dataset:
data(tikus)
tikusdat <- mvabund(tikus$abund)
year <- tikus$x[,1]

## Plot mean-variance plot for a mvabund object with a log scale (default):
meanvar.plot(tikusdat) 	

## Again but without log-transformation of axes:
meanvar.plot(tikusdat,log="") 	

## A mean-variance plot, data organised by year, 
## for 1981 and 1983 only, as in Figure 7a of Warton (2008):
is81or83 <- year==81 | year==83
meanvar.plot(tikusdat~year, subset=is81or83, col=c(1,10))