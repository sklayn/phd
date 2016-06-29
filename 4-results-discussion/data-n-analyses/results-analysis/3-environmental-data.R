### Environmental data - missing values, imputations
### (continued from community composition)

# define the working subdirectories
data.dir <- "data"
functions.dir <- "R"
save.dir <- "output"
figs.dir <- "figs"

# import necessary packages
library(ggplot2)
library(grid)
library(mice)
library(miceadds)
library(plyr)
library(reshape2)
library(vegan)
library(VIM)


## import environmental data - for the study period, as is, in order to explore
## (missing values, etc.) and decide what to do next
# env.sand <- read.csv()

## First, explore and decide how to deal with missing values

# get the proportion of missing values in each row (sample) and each column
# (variable)
apply(env.sand, 2, function(x) {sum(is.na(x))/length(x) * 100})
apply(env.sand, 1, function(x) {sum(is.na(x))/length(x) * 100})

# check the pattern of the missing data 
md.pattern(env.sand)

# another way - more visual (package VIM)
na.pattern <- aggr(env.sand, col = c("blue", "red"), numbers = TRUE, sortVars = TRUE,
                   labels = names(env.sand), cex.axis = 0.7, gap = 3)

## Here - a lot of missing values for the study period, so removing the whole rows
## would significantly decrease the size of the dataset and would probably lead to
## biased results. So, we'll try to impute the missing values from available data. 
## For the water column parameters, long-term monitoring data is available for our 
## stations, so all of it will be used. However, for the sediment parameters (and 
## some water column parameters such as O2, temperature, salinity) only short-term 
## observation data is available, so these will be imputed separately (to avoid 
## having to impute a very large set of irrelevant and possibly very improbable values). 

## Imputation of missing data (package mice, miceadds)

# import the long-term water column parameters (raw - include all available measures)  
water.param.sand <- read.csv(file.path(data.dir, "waterColumnLT_sand-imputations.csv"),
                             header = T)

# here only: rename Atia to Akin, and Maslen nos to Paraskeva - to
# have more values for the missing value imputation. These stations are located very 
# close to ours, so it is not too wrong to do this. 
levels(water.param.sand$station) <- revalue(levels(water.param.sand$station), 
                                            c("Atia"="Akin", "Maslen nos"="Paraskeva"))

# reorder factor levels for station as desired
water.param.sand$station <- reorder.factor(water.param.sand$station, 
                                           new.order = stations.sand)

# aggregate the data by year, month and station (-> average of all depths for a
# particular sampling occasion)
water.param.sand.aggr <- ddply(water.param.sand, .(year, month, station), colwise(mean, .cols = is.numeric))

# remove the depth (not relevant for these particular imputations)
water.param.sand.aggr$depth <- NULL

# check the predictor matrix. The correlation between predictor and target, and 
# the proportion of usable cases are used for its construction. 
quickpred(water.param.sand.aggr)

# impute the missing data
water.sand.imp <- mice(water.param.sand.aggr, method = "pmm", m = 100, seed = 100, printFlag = FALSE)
summary(water.sand.imp)

# extract the imputed data (long format) and subset to 2013-2014 only (the only
# ones we're interested in)
water.sand.imp.df <- complete(water.sand.imp, action = "long", include = FALSE)
water.sand.imp.subs <- subset(water.sand.imp.df, year == 2013 | year == 2014)

# sort each imputation (id contained within column .imp) by station, then years and months
water.sand.imp.subs <- ddply(water.sand.imp.subs, .(.imp), function(x) arrange(x, station, year, month))

# get rid of the (confusing) .id column
water.sand.imp.subs$.id <- NULL  


# import the sediment parameters (and other parameters for which only short-term 
# data is available)
other.env.sand <- read.csv(file.path(data.dir, "other-env_sand-imputations.csv"),
                           header = TRUE)

# reorder factor station as desired
other.env.sand$station <- reorder.factor(other.env.sand$station, 
                                         new.order = stations.sand)


# impute missing values here, too

# first, change the predictor matrix so that the total heavy metals is always the
# sum of all other heavy metals (idem for heavy metals w/o Fe).
# ini <- mice(other.env.sand, maxit = 0, printFlag = FALSE)
# meth <- ini$meth
# meth["heavy.metals.all"] <- "~I(Cu + Pb + Zn + Cd + Mn + Fe + Ni)"
# meth["heavy.metals.noFe"] <- "~I(Cu + Pb + Zn + Cd + Mn + Ni)"

# now, impute the missing values
other.env.sand.imp <- mice(other.env.sand, 
                           method = "pmm", 
                           m = 100, 
                           seed = 100, 
                           printFlag = FALSE)

summary(other.env.sand.imp)

# extract the imputed data (long format) and subset to 2013-2014 only
other.env.sand.imp.df <- complete(other.env.sand.imp, action = "long", include = FALSE)

# fix the heavy metals and heavy metals w/o Fe (doesn't seem to work before imputation)
other.env.sand.imp.subs <- ddply(other.env.sand.imp.df, 
                                 .(.imp), 
                                 transform, 
                                 heavy.metals.all = Cu + Pb + Zn + Cd + Mn + Fe + Ni, 
                                 heavy.metals.noFe = Cu + Pb + Zn + Cd + Mn + Ni)

# sort each imputation (id contained within column .imp) by station, then years and months
other.env.sand.imp.subs <- ddply(other.env.sand.imp.subs, .(.imp), function(x) arrange(x, station, year, month))

# get rid of the (confusing) .id column
other.env.sand.imp.subs$.id <- NULL  

# subset to 2013-2014 only
other.env.sand.imp.subs <- subset(other.env.sand.imp.subs, year == 2013 | year == 2014)


# combine with the water column (imputed and subsetted) data 
env.imp.all.sand <- join(water.sand.imp.subs, other.env.sand.imp.subs)


## get the long-term average data for our stations (from both imputed data sets - water column and sediment - 
## basically combine without subsetting & clean); then set aside (will perform all analyses with them, too)

# average the observations by station & get rid of no longer necessary year & month
water.imp.LT.sand <- ddply(water.sand.imp.df[, !names(water.sand.imp.df) %in% c("month", "year")], 
                           .(station), colwise(mean, .cols = is.numeric))

other.env.imp.LT.sand <- ddply(other.env.sand.imp.df[, !names(other.env.sand.imp.df) %in% c("month", "year")], 
                               .(station), colwise(mean, .cols = is.numeric))

# merge the two long-term environmental data frames
env.imp.all.LT.sand <- join(water.imp.LT.sand, other.env.imp.LT.sand, by = "station") 

# repeat the ugly hack to obtain the same number of rows for all replicates
# (copy each row 9 times so that the LT values here match the station replicates 
# in the abundance data frame, etc.)
env.imp.all.LT.sand <- env.imp.all.LT.sand[rep(seq_len(nrow(env.imp.all.LT.sand)), each = 9), ]
rownames(env.imp.all.LT.sand) <- 1:nrow(env.imp.all.LT.sand)

# save LT data and set aside for later analyses                              
saveRDS(env.imp.all.LT.sand, file.path(save.dir, "env-data-imputed-LT-clean_sand.rds"))


# clean up the no longer useful intermediate data frames, lists etc. 
rm(water.param.sand, 
   water.param.sand.aggr,
   water.sand.imp.df,
   water.sand.imp, 
   water.sand.imp.subs,
   water.imp.LT.sand,
   other.env.sand, 
   other.env.sand.imp.df,
   other.env.sand.imp, 
   other.env.sand.imp.subs, 
   other.env.imp.LT.sand)
# rm(ini, meth)

# repeat each row 3 times to match the rows of the mds object and fix row names
# (very ugly, but working, hack)
env.imp.all.sand <- env.imp.all.sand[rep(seq_len(nrow(env.imp.all.sand)), each = 3), ]
rownames(env.imp.all.sand) <- 1:nrow(env.imp.all.sand)

# SAVE THIS IF NECESSARY FOR FUTURE USE 
saveRDS(env.imp.all.sand, file.path(save.dir, "env-data-imputed-clean_sand.rds"))


########################################################################################################
### ASIDE: 
### Gradient forest analysis -> try to predict species presence by the available 
### environmental variables (package extendedForest, gradientForest). Most probably
### won't use because no time to resolve warnings, but just to keep in mind to 
### explore in the future
library(gradientForest)

# prepare the environemental data: use list of all (100) imputed values for the 
# environmental variables; combine by averaging by station 
# => equivalent to vars.for.ordisurf
env.vars.sand.gf <- ddply(env.imp.all.sand,
                          .(station, year, month), colwise(mean, .cols = is.numeric))

# repeat each row 3 times (to match the 3 replicates per station in the abundance data)
env.vars.sand.gf <- env.vars.sand.gf[rep(seq_len(nrow(env.vars.sand.gf)), each = 3), ]
rownames(env.vars.sand.gf) <- 1:nrow(env.vars.sand.gf)


# calculate the maximum number of splits for the random tree analysis
lev <- floor(log2(nrow(num.zoo.abnd.sand) * 0.368/2))

# apply the gradient forest analysis (using the subset of species present at the 
# current sites).
# Important: convert the sites x species df to matrix; exclude the first 3 columns
# (station, year and month) from the df of environmental variables; include a 
# transformation of the species abundance data (here - square root). 
gf.sand <- gradientForest(cbind(env.vars.sand.gf[-c(1:3)], as.matrix(num.zoo.abnd.sand)), 
                          predictor.vars = colnames(env.vars.sand.gf[-c(1:3)]), 
                          response.vars = colnames(num.zoo.abnd.sand), 
                          ntree = 500,
                          transform = function(x) sqrt(x),
                          compact = TRUE,
                          nbin = 201,
                          maxLevel = lev,
                          corr.threshold = 0.5
                          )

gf.sand

# plots to explore the results... 

######################################################################################################################

## PCA on environmental parameters
library(FactoMineR)
library(factoextra)

# if not installed:
# devtools::install_github("kassambara/factoextra")
# install.packages("FactoMineR", "factoextra")

# only on 2013-2014 data (with previously imputed missing values)















## Dissimilarities and environment
# PERMANOVA - multivariate ANOVA based on dissimilarities (=> adonis in vegan)

# We'll study beta diversity using the variables identified as most significant
# above (to ease computation; later if time might try others, too for kicks). 
# Beta diversity here = slope of the species-area curve (= Arrhenius z). Most 
# common interpretation: z = 0.3 - random sampling variability; higher values - 
# real systematic differences. 
beta.div.sand <- betadiver(num.zoo.abnd.sand, "z")
summary(beta.div.sand)

# check the multivariate homogeneity of groups dispersions (variances)
# => assumption of PERMANOVA
anova(betadisper(beta.div.sand, group = factors.zoo.sand$stations))

