### Environmental data - missing values, imputations
### (continued from community composition)

# define the working subdirectories
data.dir <- "data"
functions.dir <- "R"
save.dir <- "output"
figs.dir <- "figs"

# import necessary packages
library(mice)
library(miceadds)
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
water.sand.imp <- mice(water.param.sand.aggr, 
                             method = "pmm", 
                             m = 100, 
                             seed = 100, 
                             printFlag = FALSE)
summary(water.sand.imp)

# extract the imputed data (long format) 
water.sand.imp.df <- complete(water.sand.imp, action = "long", include = FALSE)

# sort each imputation (id contained within column .imp) by station, then years and months
water.sand.imp.df <- ddply(water.sand.imp.df, 
                           .(.imp), 
                           function(x) arrange(x, station, year, month))

# get rid of the (confusing) .id column
water.sand.imp.df$.id <- NULL  

## SAVE THIS TO AVOID REPEATING THE IMPUTATION PROCEDURE WHICH IS A RESOURCE HOG
saveRDS(water.sand.imp.df, file = file.path(save.dir, "water-column-imputed_clean_sand.rds"))


# import the sediment parameters (and other parameters for which only short-term 
# data is available)
other.env.sand <- read.csv(file.path(data.dir, "other-env_sand-imputations.csv"),
                           header = TRUE)

# reorder factor station as desired
other.env.sand$station <- reorder.factor(other.env.sand$station, 
                                         new.order = stations.sand)

# impute missing values here, too
other.env.sand.imp <- mice(other.env.sand, 
                           method = "pmm", 
                           m = 100, 
                           seed = 100, 
                           printFlag = FALSE)

summary(other.env.sand.imp)

# extract the imputed data (long format) and subset to 2013-2014 only
other.env.sand.imp.df <- complete(other.env.sand.imp, action = "long", include = FALSE)

# sort each imputation (id contained within column .imp) by station, then years
# and months (also unnecessary, since already sorted)
other.env.sand.imp.df <- ddply(other.env.sand.imp.df, 
                                 .(.imp), 
                                 function(x) arrange(x, station, year, month))

# get rid of the (confusing) .id column
other.env.sand.imp.df$.id <- NULL  

## SAVE THIS DATASET TOO - TO AVOID IMPUTING AGAIN
saveRDS(other.env.sand.imp.df, file = file.path(save.dir, "other-env-imputed_clean_sand.rds"))


# clean up the no longer useful intermediate data frames, lists etc. 
rm(water.param.sand, 
   water.param.sand.aggr,
   water.sand.imp,
   other.env.sand, 
   other.env.sand.imp)


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

