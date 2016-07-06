## mvabund - GLM fitting for abundance - environmental data
library(mvabund)

## turn abundance data into mvabund object = 2D matrix (easier for mvabund 
## analyses)
sand.mvabund <- mvabund(num.zoo.abnd.sand)

## visualize the data (by certain factors that we suspect might influence the 
# community structure) --> IF USING THESE, FIX PLOTS (MOVE LEGEND & ADD LEGEND TITLE, FIX MARGINS)
plot(sand.mvabund ~ as.factor(gr.dendr.sand)) # by clusters identified from the dendrogram
plot(sand.mvabund ~ env.qualit$LUSI.3000.impact)  # by level of anthropogenic pressure (LUSI)
# do it for other factors (group has to be a factor!)


## checking model assumptions
# 1. mean-variance assumption => choice of family parameter. Can be checked by 
# plotting residuals vs fits: if little pattern - the chosen mean-variance 
# assumption is plausible.
# another way: direct plotting (variance ~ mean), for each species w/n each factor
# level.

# here: using groups from dendrogram as factors. 
# NB do this plot several times: these residuals use random number generation, so 
# replicate runs will not be identical, but the pattern should stay consistent
# across replicate plots. 
plot(manyglm(sand.mvabund ~ gr.dendr.sand, family = "negative.binomial"))

meanvar.plot(sand.mvabund ~ as.factor(gr.dendr.sand), table = TRUE)
meanvar.plot(sand.mvabund, table = TRUE)


# 2. assumed relationship b/n mean abundance and environmental variables - link 
# function and formula.
# If quantitative variables included in the model -> if trend in size of residuals
# at different fitted values (e.g. U-shape,..) = violation of the log-linearity 
# assumption.
#  

## if everything looks fine, fit the model
sand.glms <- manyglm(sand.mvabund ~ gr.dendr.sand, family = "negative.binomial")

## explore the fit
plot(sand.glms)  # residuals vs fitted values
plot.manyglm(sand.glms, which = 1:4)  # all traditional lm diagnostic plots

# get specific model parameters
residuals(sand.glms)
coef(sand.glms)
fitted.values(sand.glms)

# NB summary can take a long time, depending on number of resamplings!
summary(sand.glms, test = "LR", p.uni = "adjusted", show.time = "total")


## Testing hypotheses about the community-environment association
# -> anova on manyglm result. NB Can take a long time (depending on number of resamplings!)
# NB 2 temporal autocorrelation (revisited sites!)

strata <- rep(1:3, each = 3, times = 6)  # probably no good, because also replicates at each site
  
sand.anova.glm <- anova.manyglm(sand.glms, 
                                test = "LR", 
                                block = strata,
                                nBoot = 100,
                                p.uni = "adjusted", 
                                show.time = "total")

# Relative taxon contribution to patterns: LR statistic - a measure of strength 
# of individual taxon contributions

# retrieve and plot species with the highest contribution to the patterns tested -> DOESN'T WORK IN THIS 
uniSorted <- sort(sand.anova.glm$uni.test, index.return = TRUE, decreasing = TRUE)
plot(sand.mvabund[ , uniSorted$ix[1:10]]) 



# advantage of multivariate analysis - greater power to detect patterns when 
# analysing all species simultaneously than when looking for a pattern specifically
# in each species.

