## Environmental data 

library(mice)
library(miceadds)

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

# impute the missing data
env.imp <- mice(env.sand, method = "pmm", m = 5)
summary(env.imp)


# complete the observed data using the first (or any other you want) of the 
# imputed datasets -> BETTER DO ANALYSIS ON ALL IMPUTED DATASETS THEN POOL THE
# FINAL RESULTS!
env.completed <- complete(env.imp, 1)

# compare the distributions of the original and the imputed data 
stripplot(env.imp)

#####################################################################################################



# get the long-term water column parameters (raw - include all available measures)  
water.param <- read.csv(file.path(data.dir, "waterColumnLT-imputations.csv"),
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
other.env <- read.csv(file.path(data.dir, "other-env-imputations.csv"),
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

# (here only, for testing) fix the total heavy metals and heavy metals w/o Fe
# (sums of the other heavy metals, and if one is missing - all are missing)
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
 
# sort each data frame in the list by p-value
env.mds.sign <- lapply(env.mds.sign, function(x) arrange(x, pvals, desc(r2)))

# count the number of times a variable occurs at a particular position   
# only get vars + row names (variable & position) from each table
x <- lapply(env.mds.sign, function(x) {
                            # add a variable for the position
                            x$pos <- row.names(x)
                            # subset to get desired variables only
                            x <- subset(x, select = c(vars, pos))
                            return(x)
                          })

# convert list of data frames to 1 big data frame
x <- do.call("rbind", x)

# count the total number of times a variable appears as significant in the 
# imputed datasets
y1 <- arrange(count(x, vars = "vars"), desc(freq))

# convert pos to numeric to facilitate sorting & sort data frame on position 
# and (descending) frequency
y$pos <- as.numeric(y$pos)

# count the occurrences of each variable in each position
y <- count(x, vars = c("vars", "pos"))
y <- arrange(y, pos, desc(freq))
y

# plot each variable's frequency by position (variables as labels)
with(y, text(pos, freq, labels = vars))


## use the most significant environmental variables to plot as surfaces
## (ordisurf in package vegan) -- MAKE INTO A FUNCTION!
ordi.tst <- ordisurf(mds.sand ~ num.env.imp.all[[1]]$Secchi.depth)

# get data for plotting
ordi.grid <- ordi.tst$grid #extracts the ordisurf object
str(ordi.grid) #it's a list though - cannot be plotted as is
ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi.mite$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.mite.na <- data.frame(na.omit(ordi.mite)) #gets rid of the nas
ordi.mite.na #looks ready for plotting!


#data for plotting
##NMDS points
NMDS.data <- env.imp.all[[1]] #there are other ways of doing this. But this is the way I do it for ease of plotting
NMDS.data$NMDS1 <- mds.sand$points[ ,1] #this puts the NMDS scores for the plots into a new dataframe. you could put them into an existing one if you preferred.
NMDS.data$NMDS2 <- mds.sand$points[ ,2]


## Plotting in ggplot2
ordisurf.ggplot <- ggplot(NMDS.data, aes(x = NMDS1, y = NMDS2)) +
                      geom_contour(data = ordi.mite.na, aes(x = x, y = y, z = z, colour = ..level..),
                          bins = 8) + #can change the binwidth depending on how many contours you want
                      geom_point(size = 3, alpha = 0.8) + #plots the NMDS points
                      theme_bw() +
                      labs(colour = "Secchi depth") + #another way to set the labels, in this case, for the colour legend
                      scale_colour_gradient(high = "darkblue", low = "skyblue") + #here we set the high and low of the colour scale.  Can delete to go back to the standard blue, or specify others
                      theme(legend.key = element_blank(),  #removes the box around each legend item
                            legend.position = "bottom", #legend at the bottom
                            legend.direction = "horizontal",
                            legend.box = "horizontal",
                            legend.box.just = "centre")
ordisurf.ggplot