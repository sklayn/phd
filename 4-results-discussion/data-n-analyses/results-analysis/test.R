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

# check if transforming the data helps with the mean-varaiance relationship - 
# no, not at all, actually makes it worse. # This happens because of all the 
# zeros in the data, and is typical of this type of data.
meanvar.plot(log(sand.mvabund + 1) ~ as.factor(gr.dendr.sand))
meanvar.plot(sqrt(sand.mvabund) ~ as.factor(gr.dendr.sand))


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

repl.blocks <- rep(1:3, each = 3, times = 6)  # probably no good, because also replicates at each site
  
sand.anova.glm <- anova.manyglm(sand.glms, 
                                test = "LR",
                                p.uni = "adjusted", 
                                show.time = "all")

print(sand.anova.glm$table)


## Relative taxon contribution to patterns: LR statistic - a measure of strength 
# of individual taxon contributions. 
# (LR expresses how many times more likely the data are under one model than the 
# other. This likelihood ratio, or equivalently its logarithm, can then be used 
# to compute a p-value, or compared to a critical value to decide whether to 
# reject the null model in favour of the alternative model.)


## retrieve and plot (top 10) species with the highest contribution to the patterns tested

# get 2nd row from table of univariate test stats (= change in deviance due to groups)
uniSorted <- sort(sand.anova.glm$uni.test[2, ], decreasing = TRUE, index.return = TRUE) 

# get those species' names and plot them
dimnames(sand.mvabund)[[2]][uniSorted$ix[1:10]]
plot(sand.mvabund ~ as.factor(gr.dendr.sand), var.subset = uniSorted$ix[1:10]) 

# get the percentage of change in deviance due to this "top ten" => only about 
# 50%, so cannot focus on just them if we want the whole story
sum(sand.anova.glm$uni.test[2, uniSorted$ix[1:10]]) / sum(sand.anova.glm$uni.test[2, ])

# with 50 species, most of the change in deviance is explained
sum(sand.anova.glm$uni.test[2, uniSorted$ix[1:50]]) / sum(sand.anova.glm$uni.test[2, ])


# advantage of multivariate analysis - greater power to detect patterns when 
# analysing all species simultaneously than when looking for a pattern specifically
# in each species.



## Relate the species abundances to environmental data

# prepare the environmental data (=> reduced nb of variables, because otherwise
# fries all my computers) 
env.all.sand.glm <- subset(env.all.sand.envfit, select = pruned.cors.vars)

# get rid of singletons - rarely informative AND might reduce computation time
sand.mvabund.reduced <- num.zoo.abnd.sand[, which(specnumber(num.zoo.abnd.sand, MARGIN = 2) > 1)]
sand.mvabund.reduced <- as.mvabund(sand.mvabund.reduced)

# perform the analysis
sp.env.glm.sand <- manyglm(paste("sand.mvabund.reduced ~", paste(names(env.all.sand.glm), collapse = "+")),
                           data = env.all.sand.glm, 
                           family = "negative.binomial")

anova.sp.env.glm.sand <- anova.manyglm(sp.env.glm.sand, 
                                       test = "LR",
                                       p.uni = "adjusted",
                                       nBoot = 10,
                                       show.time = "all" 
                                       )
print(anova.sp.env.glm.sand$table)




## see which species respond differently to different environmental parameters
## (= fit single predictive model for all species at all sites, but w/o attempting
## to explain the different responses using traits - the species ID is used in 
## place of a traits matrix). 
## Use a LASSO penalty (by setting the term method) when applying, to simplify 
## interpretation: does automatic model selection, setting to zero any interaction
## coefficients that don’t help reduce BIC.

## try with different factors (should be factors!): groups in the MDS / LUSI / sth else..


## NB should use reduced env.variables dataset => for ex. important variables 
## identified through PCA, CCA, envfit,.. -> OR RUN AT HOME W/ EVERYTHING AND THEN COMPARE
ft.sp.env <- traitglm(sand.mvabund.reduced, 
                      env.all.sand.glm, 
                      method = "glm1path")

ft.sp.env$fourth

# plot this -> make 2 subplots w/ ~50 species each, or will be unreadable!
a <- max(abs(ft.sp.env$fourth.corner))
colort <- colorRampPalette(c("blue","white","red")) 
plot.spp <- levelplot(t(as.matrix(ft.sp.env$fourth.corner)), xlab = "Environmental Variables",
                     ylab = "Species", col.regions = colort(100), at = seq(-a, a, length = 100),
                     scales = list(x = list(rot = 45)))
print(plot.spp)




###################################################################################################################
###################################################################################################################
### TESTINGGG

## average abundances by replicate (=> 18 rows = 6 stations x 3 samplings)
tst.mvabund <- cbind(num.zoo.abnd.sand, factors.zoo.sand)
tst.mvabund$replicates <-  rep(1:3, each = 3, times = 6)
tst.mvabund <- ddply(tst.mvabund, .(replicates, stations), colwise(mean, .cols = is.numeric))
tst.mvabund <- arrange(tst.mvabund, stations)
head(tst.mvabund)

tst.sand.mvabund <- mvabund(tst.mvabund[-c(1:2)])

# make new grouping factors ot correspond to order of rows in the new matrix:
# dendr. groups
tst.groups <- c(rep(1, 6), rep(2, 3), rep(3, 3), rep(4, 3), rep(3, 3))

# plot
plot(tst.sand.mvabund ~ as.factor(tst.groups)) 

## check model assumptions
plot(manyglm(tst.sand.mvabund ~ tst.groups, family = "negative.binomial"))
meanvar.plot(tst.sand.mvabund ~ as.factor(tst.groups), table = TRUE)

## look ok; fit model
tst.sand.glms <- manyglm(tst.sand.mvabund ~ tst.groups, family = "negative.binomial")

## explore the fit
plot(tst.sand.glms)  # residuals vs fitted values
plot.manyglm(tst.sand.glms, which = 1:4)  # all traditional lm diagnostic plots

residuals(tst.sand.glms)
coef(tst.sand.glms)
fitted.values(tst.sand.glms)

summary(tst.sand.glms, test = "LR", nBoot = 30, p.uni = "adjusted", show.time = "total")


## Testing hypotheses about the community-environment association
# -> anova on manyglm result. NB Can take a long time (depending on number of resamplings!)
# NB 2 temporal autocorrelation (revisited sites!)

tst.strata <- rep(1:3, times = 6)  

tst.sand.anova.glm <- anova.manyglm(tst.sand.glms, 
                                    test = "LR", 
                                    nBoot = 30,
                                    p.uni = "adjusted", 
                                    show.time = "total")

# retrieve and plot species with the highest contribution to the patterns tested -> DOESN'T WORK IN THIS 
tst.uniSorted <- sort(tst.sand.anova.glm$uni.test, index.return = TRUE, decreasing = TRUE)
plot(tst.sand.mvabund[ , tst.uniSorted$ix[1:10]]) 



######################################################################################################################