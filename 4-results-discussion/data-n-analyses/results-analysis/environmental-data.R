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

# aggregate the data by year, month and station (-> average of all depths for a
# particular sampling occasion)
water.param.aggr <- ddply(water.param, .(year, month, station), colwise(mean, .cols = is.numeric))

# remove the depth (not relevant for these particular imputations)
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

# sort by station, then years and months
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

# clean up the no longer useful intermediate data frames, lists etc. 
rm(water.param, water.param.aggr, water.imp, water.imp.subs, other.env, other.env.imp, other.env.imp.subs)


# (here only) fix the total heavy metals and heavy metals w/o Fe - these variables 
# don't need to be imputed, because they are sums of the other individual heavy metals, 
# and if one of those is missing, all are missing (= all are imputed).
env.imp.all <- lapply(env.imp.all, function(x) { 
                                      x["heavy.metals.all"] <- x$Cu + x$Pb + x$Zn + x$Cd + x$Mn + x$Fe + x$Ni
                                      x["heavy.metals.noFe"] <- x$Cu + x$Pb + x$Zn + x$Cd + x$Mn + x$Ni
                                      return(x)
                                   })

# order data frames by station and date
env.imp.all <- lapply(env.imp.all, function(x) arrange(x, station, year, month))

# repeat each row 3 times to match the rows of the mds object and fix row names
# (very ugly hack)
env.imp.all <- lapply(env.imp.all, function(x) x[rep(seq_len(nrow(x)), each=3),])
env.imp.all <- lapply(env.imp.all, function(x) {rownames(x) <- 1:nrow(x); x})


# run envfit (vegan) on the imputed datasets (env. variables subset only)
env.mds.sand <- lapply(env.imp.all, function(x) {
                                            # subset each data frame to only the variable columns
                                            envfit(mds.sand, 
                                                   subset(x, select = -c(station, month, year)), 
                                                   permutations = 999)
                                        })

sign_vars_pos_envfit <- function(envfit.obj) {
  ## helper function to extract the most significant (p < 0.05) variables
  ## and their position (based on p and r2).
  
  # extract the most significant variables (p < 0.05) from the fit
  envfit.sign.vars <- extract_envfit_scores(envfit.obj, p = 0.05, r2 = T)
  
  # sort the resulting data frame by p-value and descending r-squared
  envfit.sign.vars <- arrange(envfit.sign.vars, pvals, desc(r2))
  
  # count the number of times a variable occurs at a particular position:   
  # only get vars + row names (variable & position) from the data frame
  envfit.sign.vars$pos <- row.names(envfit.sign.vars)
  sign.vars.pos <- subset(envfit.sign.vars, select = c(vars, pos))
  
  # covert pos to numeric for easier sorting
  sign.vars.pos$pos <- as.numeric(sign.vars.pos$pos)
  
  return(sign.vars.pos)  
} 

# extract the most significant variables and positions from the list of envfits
env.sign.count <- lapply(env.mds.sand, sign_vars_pos_envfit)

sign_vars_freq <- function(env.vars.count, target.freq) {
  ## completely useless helper function, to be applied only here, for examining
  ## the frequency of occurence of environemental variables across a list of 
  ## variables and their positions obtained from the respective envfits. 
  ## Extracts the names of the final variables - to be plotted as surfaces or
  ## vectors over an mds.
  ## Input is a list of data frames (each - columns vars and pos).
  
  # convert list of data frames to 1 big data frame
  sign.count.df <- do.call("rbind", env.vars.count)
  
  # 1. count the total number of times a variable appears as significant across the 
  # imputed datasets
  sign.count.all <- arrange(count(sign.count.df, vars = "vars"), desc(freq))

  # 2. count the occurrences of each variable in each position
  sign.count.pos <- count(sign.count.df, vars = c("vars", "pos"))
  sign.count.pos <- arrange(sign.count.pos, pos, desc(freq))

  # examine the distribution of the variables and select the most significant 
  # (the ones that appear the most often, i.e. with a frequency >= target frequency)
  # at the topmost positions after the envfit procedure)
  sign.vars <- unique(subset(sign.count.pos, freq >= target.freq)$vars)
  
  # drop the unused factor levels & reorder the factor (order is important here = 
  # decreasing significance)
  sign.vars <- droplevels(sign.vars)
  sign.vars <- reorder.factor(sign.vars, new.order = unique(sign.vars))
  
  # return the variables as well as all counts in a list, just in case
  sign.vars.list <- list("sign.vars" = sign.vars, 
                         "sign.freq" = sign.count.all, 
                         "sign.pos.freq" = sign.count.pos) 
  return(sign.vars.list)
}

# get the most significant environmental variables according to the envfits 
# performed on the imputed datasets
sign.vars <- sign_vars_freq(env.sign.count, target.freq = 10)

# remove workspace clutter 
rm(env.sign.count)

## PLOTTING VARIABLES AS SURFACES
# get the original numeric values of the environmental variables chosen for plotting
# from all the imputed datasets

# only use env. variables from each df
vars.for.ordisurf <- lapply(lapply(env.imp.all, function(y) subset(y, select = -c(station, month, year))), # only env.variables from each df
                            function(x) { 
                              x <- subset(x, select = names(x) %in% levels(sign.vars[[1]])) # subset according to vector of chosen sign. variables
                              # add a numeric identifier to be used later for aggregating
                              # to avoid having to reorder later
                              x$id <- rownames(x)
                              x$id <- as.numeric(x$id)
                              return(x)
                            }) 


# combine into a single data frame and average each variable by id (= replicate)
vars.for.ordisurf <- do.call("rbind", vars.for.ordisurf) 
vars.for.ordisurf <- ddply(vars.for.ordisurf, .(id), colwise(mean)) 
 
# apply ordisurf sequentially to all environmental variables (get back a list of
# ordisurf objects where each element is an environmental variable)
ordisurf.list.all <- apply(vars.for.ordisurf, 
                       MARGIN = 2, 
                       FUN = function(x) ordi <- ordisurf(mds.sand ~ x, plot = F)) 

# get rid of the first element (id), which serves no purpose
ordisurf.list.all$id <- NULL

# check out the summaries of the fits
lapply(ordisurf.list.all, summary)

# drop O2.average - mostly the same as O2.bottom (make new list, in case it's 
# actually needed later)
ordisurf.list.noO2 <- ordisurf.list.all
ordisurf.list.noO2$O2.average <- NULL

# rearrange list to have the plots in the desired order.
# Here: on the first row of the plot will be the sediment parameters, on the 
# second - the water column parameters, and on the third - the heavy metals.

ordisurf.list.noO2 <- ordisurf.list.noO2[c("sorting", "mean.grain.size", "org.matter", "silt.clay", 
                                   "O2.bottom", "Secchi.depth", "seston", "salinity",
                                   "dist.innermost", "depth", "heavy.metals.noFe", "Pb",
                                   "Ni", "Cu")]

var.labels <- c("sorting", "mean grain size", "organic matter", "silt-clay", 
                 "O2 bottom", "Secchi depth", "seston", "salinity",  
                 "distance to innermost station", "depth", "heavy metals (no Fe)", "Pb",
                 "Ni", "Cu")

# set the file name and properties for the output graph
pdf(file = file.path(figs.dir, "mds_ordisurf_sand_most_sign_vars.pdf"), 
    paper = "a4",
    width = 12,
    height = 12,
    useDingbats = F)

# modify par to fit all plots on one page (here, 4 plots per row)
par(mfrow = c(5, 3))

# plot all variables, using the custom plot_mds_ordisurf function, and adding the
# corresponding main title (variable name) on each subplot
mapply(function(m, n) {
          plot_mds_ordisurf(mds.sand, m)
          title(main = n, col.main = "grey28")
       }, 
       ordisurf.list.noO2, 
       var.labels
       )

dev.off()

# return the graphics device to the original settings
par(mfrow = c(1, 1))

# clean up workspace 
rm(ordisurf.list.noO2, ordisurf.list.all, var.labels, vars.for.ordisurf)


### FIX FROM HERE
### Gradient forest analysis -> try to predict species presence by the available 
### environmental variables (package extendedForest, gradientForest)
library(gradientForest)

# prepare the environemental data: use list of all (100) imputed values for the 
# environmental variables; combine by averaging by station
num.env.imp.all <- lapply(lapply(env.imp.all, function(y) subset(y, select = -c(station, month, year))), # only env.variables from each df 
                                  function(x) { 
                                  # add a numeric identifier to be used later for aggregating
                                  # to avoid having to reorder later
                                  x$id <- rownames(x)
                                  x$id <- as.numeric(x$id)
                                  return(x)
                                  })
env.vars.gradientForest <- do.call("rbind", num.env.imp.all)
env.vars.gradientForest <- ddply(env.vars.gradientForest, .(id), colwise(mean))

# calculate the maximum number of splits for the random tree analysis
lev <- floor(log2(nrow(community.sand) * 0.368/2))

# apply the gradient forest analysis (using the subset of species present at the 
# current sites).
# Important: convert the sites x species df to matrix; exclude the id column from
# the df of environmental variables; include a transformation of the species abundance
# data (here - square root). 
gf.sand <- gradientForest(cbind(num.values.for.ordiplot.mean[-1], as.matrix(community.sand)), 
                          predictor.vars = colnames(num.values.for.ordiplot.mean[-1]), 
                          response.vars = colnames(community.sand), 
                          ntree = 500,
                          transform = function(x) sqrt(x),
                          compact = T,
                          nbin = 201,
                          maxLevel = lev,
                          corr.threshold = 0.5
                          )
gf.sand
