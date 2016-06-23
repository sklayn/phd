### Environmental data & ordination
### (continued from community composition)

library(ggplot2)
library(grid)
library(mice)
library(miceadds)
library(plyr)
library(reshape2)
library(vegan)

## import environmental data


## First, explore and deal with missing values

# get the proportion of missing values in each row (sample) and each column
# (variable)
apply(env.sand, 2, function(x) {sum(is.na(x))/length(x) * 100})
apply(env.sand, 1, function(x) {sum(is.na(x))/length(x) * 100})

# check the pattern of the missing data 
md.pattern(env.sand)

# another way - more visual
library(VIM)
na.pattern <- aggr(env.sand, col = c("blue", "red"), numbers = T, sortVars = T,
                   labels = names(env.sand), cex.axis = 0.7, gap = 3)


#####################################################################################################



# get the long-term water column parameters (raw - include all available measures)  
water.param <- read.csv(file.path(data.dir, "waterColumnLT_sand-imputations.csv"),
                        header = T)

# here only, for testing: rename Atia to Akin, and Maslen nos to Paraskeva - to
# have more values for the missing value imputation
levels(water.param$station) <- revalue(levels(water.param$station), 
                                       c("Atia"="Akin", "Maslen nos"="Paraskeva"))

# reorder factor levels for station as desired
water.param$station <- reorder.factor(water.param$station, 
                                      new.order = stations.sand)

# aggregate the data by month and station (average of all depths)
water.param.aggr <- ddply(water.param, .(year, month, station), colwise(mean, .cols = is.numeric))

# remove the depth (not relevant for these imputations)
water.param.aggr$depth <- NULL

# check the predictor matrix. The correlation between predictor and target, and 
# the proportion of usable cases are used for its construction. 
quickpred(water.param.aggr)

# impute the missing data
water.imp <- mice(water.param.aggr, method = "pmm", m = 100, seed = 100)
summary(water.imp)

# subset the imputed data to only 2013-2014 (the only ones we're interested in)
water.imp.subs <- subset_datlist(water.imp, index = 1:100, 
                                 subset = water.imp[[2]]$year == 2013 | water.imp[[2]]$year == 2014)

# sort by station names
water.imp.subs <- lapply(water.imp.subs, function(x) arrange(x, station, year, month))


# import the sediment parameters (and other parameters for which only short-term 
# data is available)
other.env <- read.csv(file.path(data.dir, "other-env_sand-imputations.csv"),
                      header = T)

# reorder factor station as desired
other.env$station <- reorder.factor(other.env$station, 
                                    new.order = stations.sand)


# impute missing values here, too
other.env.imp <- mice(other.env, method = "pmm", m = 100, seed = 100)
summary(other.env.imp)

other.env.imp.subs <- subset_datlist(other.env.imp, index = 1:100, 
                                     subset = other.env.imp[[2]]$year == 2013 | other.env.imp[[2]]$year == 2014)

# combine with the water column (imputed and subsetted) data 
env.imp.all <- mapply(function(x, y) merge(x, y, by = c("station", "month", "year")), 
                      water.imp.subs, 
                      other.env.imp.subs, 
                      SIMPLIFY = F)

# (here only) fix the total heavy metals and heavy metals w/o Fe - these variables 
# don't need to be imputed, because they are sums of the other individual heavy metals, 
# and if one of those is missing, all are missing (= all are imputed).
env.imp.all <- lapply(env.imp.all, function(x) { 
                                      x["heavy.metals.all"] <- x$Cu + x$Pb + x$Zn + x$Cd + x$Mn + x$Fe + x$Ni
                                      x["heavy.metals.noFe"] <- x$Cu + x$Pb + x$Zn + x$Cd + x$Mn + x$Ni
                                      x 
                                   })

# order data frames by station
env.imp.all <- lapply(env.imp.all, function(x) arrange(x, station, year, month))

# repeat each row 3 times to match the rows of the mds object and fix row names
# (very ugly hack)
env.imp.all <- lapply(env.imp.all, function(x) x[rep(seq_len(nrow(x)), each=3),])
env.imp.all <- lapply(env.imp.all, function(x) {rownames(x) <- 1:nrow(x); x})


# subset to (numeric) variables only - get rid of the factors in the first 3 columns
num.env.imp.all <- lapply(env.imp.all, function(x) {x <- subset(x, select = -c(station, month, year)); x})

# run envfit (vegan) on the imputed datasets
env.mds.sand <- lapply(num.env.imp.all, function(x) envfit(mds.sand, x, permutations = 999))

# extract the most significant variables (p < 0.05) from each fit

extract_envfit_scores <- function(envfit.obj, pval = 0.05, r2 = FALSE) {
  # helper function to extract NMDS vector scores from envfit objects according 
  # to the specified p-value (p < 0.05 by default). Optionally also extracts 
  # the r2 value.
 
  # extract the NMDS scores for the fitted environmental parameters 
  vector.scrs <- as.data.frame(scores(envfit.obj, display = "vectors"))
  
  # convert all the variables characterizing the vectors in the envfit object 
  # to a list 
  list.scrs <- as.list(envfit.obj$vectors)

  if(r2) {
    # get r2 values 
    vector.scrs$r2 <- list.scrs$r    
  }
  
  # get the p-values and subset the data frame according to them
  vector.scrs$pvals <- list.scrs$pvals
  sign.vectors <- subset(vector.scrs, pvals < pval)
  
  # clean up a little: get rid of the row names - add them as a variable in their
  # own column instead
  sign.vectors <- cbind(vars = rownames(sign.vectors), sign.vectors)
  rownames(sign.vectors) <- 1:nrow(sign.vectors)
  
  return(sign.vectors)
}

env.mds.sign <- lapply(env.mds.sand, extract_envfit_scores, r2 = T)
 
# sort each data frame in the list by p-value and descending r-squared
env.mds.sign <- lapply(env.mds.sign, function(x) arrange(x, pvals, desc(r2)))

# count the number of times a variable occurs at a particular position:   
# only get vars + row names (variable & position) from each table
env.sign.count <- lapply(env.mds.sign, function(x) {
                            # add a variable for the position
                            x$pos <- row.names(x)
                            # subset to get desired variables only
                            x <- subset(x, select = c(vars, pos))
                            return(x)
                          })

# convert list of data frames to 1 big data frame
env.sign.count.df <- do.call("rbind", env.sign.count)

# convert pos to numeric to be able to sort data frame on position and
# (descending) frequency
env.sign.count.df$pos <- as.numeric(env.sign.count.df$pos)

# 1. count the total number of times a variable appears as significant across the 
# imputed datasets
env.sign.count.all <- arrange(count(env.sign.count.df, vars = "vars"), desc(freq))

# 2. count the occurrences of each variable in each position
env.sign.count.pos <- count(env.sign.count.df, vars = c("vars", "pos"))
env.sign.count.pos <- arrange(env.sign.count.pos, pos, desc(freq))

# 3. get the data for plotting the environmental variables

# examine the distribution of the variables and select the most significant (here - 
# the ones that appear the most often (here - with a frequency >= 10) at the topmost 
# positions after the envfit procedure)
env.sign.count.pos
sign.vars <- unique(subset(env.sign.count.pos, freq >= 10)$vars)

# drop the unused factor levels & reorder the factor (order is important here = 
# decreasing significance)
sign.vars <- droplevels(sign.vars)
sign.vars <- reorder.factor(sign.vars, new.order = unique(sign.vars))


## PLOTTING VARIABLEs AS VECTORS
# go back to the original list of imputed data frames to get the coordinates of these 
# variables (NMDS1 and NMDS2) for plotting over the zoo data ordination with ordisurf 
# (vegan) or as vectors -> average over all datasets
env.for.ordiplot <- do.call("rbind", env.mds.sign)

# subset to only the previously identified variables of interest
env.for.ordiplot <- subset(env.for.ordiplot, vars %in% sign.vars)
env.for.ordiplot$vars <- droplevels(env.for.ordiplot$vars)
env.for.ordiplot$vars <- reorder.factor(env.for.ordiplot$vars, new.order = sign.vars)

# average the coordinates (NMDS1 and NMDS2) for each parameter - ARE WE ALLOWED TO DO THAT? SEE WHAT HAPPENS..
env.for.ordiplot.mean <- ddply(env.for.ordiplot, 
                               .(vars), 
                               summarize, 
                               NMDS1.m = mean(NMDS1), 
                               NMDS2.m = mean(NMDS2))

# plot as envfit vectors


## PLOTTING VARIABLES AS SURFACES
# get the original (imputed) numeric values of the variables chosen for plotting
# (repeat subsetting procedure, this time selecting columns = variables)
num.values.for.ordiplot <- lapply(num.env.imp.all, 
                                  function(x) { x <- subset(x, select = names(x) %in% sign.vars)
                                                # add a numeric identifier to be used later for aggregating
                                                # to avoid having to reorder later
                                                x$id <- rownames(x)
                                                x$id <- as.numeric(x$id)
                                                return(x)
                                  })

# combine into a single data frame, averaging each variable by id (= replicate)
num.values.for.ordiplot <- do.call("rbind", num.values.for.ordiplot)
num.values.for.ordiplot.mean <- ddply(num.values.for.ordiplot, .(id), colwise(mean))


# apply ordisurf sequentially to all environmental variables (get back a list of
# ordisurf objects where each element is an environmental variable)
ordi.list.all <- apply(num.values.for.ordiplot.mean, 
                       MARGIN = 2, 
                       FUN = function(x) ordi <- ordisurf(mds.sand ~ x, plot = F))

# get rid of the first element (id), which serves no purpose
ordi.list.all$id <- NULL

# check out the summaries of the fits
lapply(ordi.list.all, summary)


plot_mds_ordisurf <- function(mds.obj, 
                              ordisurf.obj, 
                              st.labels = NULL,
                              col.isolines = c("blue", "red"), 
                              col.points = "grey28",
                              col.others = "grey28"
                              ) {
  ## plot the mds and the ordisurf objects for all input variables in base 
  ## graphics, because the isolines with labels are better-looking than in 
  ## ggplot2 or lattice.
  ## 

  # create an empty plot, without axes or axis labels
  plot(mds.obj, display = "sites", type = "n", 
       axes = F,
       ann = F
       )
  
  if(is.null(st.labels)) {
    # display the stations as points, if labels are not provided as input
    points(mds.obj, display = "sites", col = col.points, cex = 0.7)  

  } else {
    # display the stations as text, using the labels provided 
    text(mds.obj, display = "sites", labels = st.labels, col = col.points, cex = 0.7)
  }
  
  # create custom color palette for plotting the isolines
  isoline.cols <- colorRampPalette(col.isolines)
  
  # get the (max) number of isolines that will be plotted automatically, to make
  # sure the color ramp covers the whole extent of the isolines. Function pretty()
  # is called under the hood by the plot.ordisurf methods for determining the 
  # best number and level of those lines.
  ncols <- length(pretty(ordisurf.obj$grid$z, n = 10))
  
  # overlay the ordisurf object 
  plot(ordisurf.obj, add = T, nlevels = 10, col = isoline.cols(ncols))
  
  # add in all the extra plot elements
  axis(1, col.axis = col.others, col = col.others, col.ticks = NULL)
  axis(2, col.axis = col.others, col = col.others, col.ticks = NULL)
  box(col = col.others)
  title(xlab = "NMDS1", ylab = "NMDS2", col.lab = col.others)

}


# drop O2.average - mostly the same as O2.bottom (make new list, in case it's 
# actually needed later)
ordi.list.noO2 <- ordi.list.all
ordi.list.noO2$O2.average <- NULL

# rearrange list to have the plots in the desired order.
# Here: on the first row of the plot will be the sediment parameters, on the 
# second - the water column parameters, and on the third - the heavy metals.
# In each row - variables left to right according to decreasing frequency in the 
# list of imputed datasets.
ordi.list.noO2 <- ordi.list.noO2[c("sorting", "mean.grain.size", "org.matter", "silt.clay", 
                                   "Secchi.depth", "O2.bottom", "chl.a", "salinity", 
                                   "heavy.metals.noFe", "Pb", "Mn", "Ni")]

var.labels <- c("sorting", "mean grain size", "organic matter", "silt-clay", 
                 "Secchi depth", "O2 bottom", "chl-a", "salinity", 
                 "heavy metals (no Fe)", "Pb", "Mn", "Ni")

# set the file name and properties for the output graph
pdf(file = file.path(figs.dir, "mds_ordisurf_sand_most_sign_vars.pdf"), 
    paper = "a4r",
    width = 12,
    height = 12,
    useDingbats = F)

# modify par to fit all plots on one page (here, 4 plots per row)
par(mfrow = c(3, 4))

# plot all variables, using the new function and adding the main title on each 
# subplot
mapply(function(m, n) {
          plot_mds_ordisurf(mds.sand, m)
          title(main = n, col.main = "grey28")
       }, 
       ordi.list.noO2, 
       var.labels
       )

dev.off()

# return the graphics device to the original settings
par(mfrow = c(1, 1))


