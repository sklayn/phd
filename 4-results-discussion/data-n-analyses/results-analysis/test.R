## mvabund - GLM fitting for abundance - environmental data
library(mvabund)

## turn abundance data into mvabund object = 2D matrix (easier for mvabund 
## analyses)
sand.mvabund <- mvabund(num.zoo.abnd.sand)

## turn dendrogram groups into a factor
gr.dendr.sand <- as.factor(gr.dendr.sand) 

## visualize the data (by certain factors that we suspect might influence the 
# community structure) --> IF USING THESE, FIX PLOTS (MOVE LEGEND & ADD LEGEND TITLE, FIX MARGINS)
## NB CHECK THAT YOU HAVE FACTORS AS GROUPS!
class(gr.dendr.sand)
class(env.qualit$LUSI.3000.impact)

plot(sand.mvabund ~ gr.dendr.sand) # by clusters identified from the dendrogram
plot(sand.mvabund ~ env.qualit$LUSI.3000.impact)  # by level of anthropogenic pressure (LUSI)
# do it for other factors (group has to be a factor!)



## checking model assumptions
# 1. mean-variance assumption => choice of family parameter. Can be checked by 
# plotting residuals vs fits: if little pattern - the chosen mean-variance 
# assumption is plausible.
# another way: direct plotting (variance ~ mean), for each species w/n each factor
# level.

# 1. Groups from MDS/dendrogram as grouping factor 
# NB do this plot several times: these residuals use random number generation, so 
# replicate runs will not be identical, but the pattern should stay consistent
# across replicate plots.
class(gr.dendr.sand)  # first make sure your factor is a factor! 

plot(manyglm(sand.mvabund ~ gr.dendr.sand, family = "negative.binomial"))

meanvar.plot(sand.mvabund ~ gr.dendr.sand, table = TRUE)
meanvar.plot(sand.mvabund, table = TRUE)

# check if transforming the data helps with the mean-varaiance relationship - 
# no, not at all, actually makes it worse. This happens because of all the 
# zeros in the data, and is typical of this type of data.
meanvar.plot(log(sand.mvabund + 1) ~ gr.dendr.sand)
meanvar.plot(sqrt(sand.mvabund) ~ gr.dendr.sand)


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
coef(sand.glms)
residuals(sand.glms)
fitted.values(sand.glms)

# NB summary can take a long time, depending on number of resamplings!
sand.glms.summ <- summary(sand.glms, test = "LR", 
                          p.uni = "adjusted", 
                          nBoot = 50, # limit the number of permutations for now, otherwise computer will choke
                          show.time = "all")


## Testing hypotheses about the community-environment association
# -> anova on manyglm result. 
# NB can take a long time (depending on number of resamplings!)
# NB2 no iid for rows - revisited sites!

# using site IDs as blocks, to ensure valid permutations
tst.sites <- rep(1:6, each = 9)

sand.anova.glm <- anova.manyglm(sand.glms, 
                                test = "LR",
                                nBoot = 50, # limited number of permutations for now  
                                block = tst.sites,
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
plot(sand.mvabund ~ gr.dendr.sand, var.subset = uniSorted$ix[1:10]) 

# get the percentage of change in deviance due to this "top ten" => only about 
# 28%, so cannot focus on just them if we want the whole story
sum(sand.anova.glm$uni.test[2, uniSorted$ix[1:10]]) / sum(sand.anova.glm$uni.test[2, ])

# with 50 species, most of the change in deviance is explained (76%)
sum(sand.anova.glm$uni.test[2, uniSorted$ix[1:50]]) / sum(sand.anova.glm$uni.test[2, ])


# advantage of multivariate analysis - greater power to detect patterns when 
# analysing all species simultaneously than when looking for a pattern specifically
# in each species.



### Relate the species abundances to environmental data  #####

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
                                       show.time = "all")

print(anova.sp.env.glm.sand$table)

## try with different factors (should be factors!): groups in the MDS / LUSI / sth else..
# LUSI levels = measure of anthropogenic impact from coastal sources
sand.glm.lusi <- manyglm(sand.mvabund.reduced ~ env.qualit$LUSI.3000.impact, 
                         family = "negative.binomial")

# check the fit (several times, because randomness)
# .... looks fine
plot(sand.glm.lusi)

sand.glm.lusi.summ <- summary(sand.glm.lusi, 
                              test = "LR", 
                              p.uni = "adjusted",
                              nBoot = 20, 
                              show.time = "all")

# check the coefficients
coef(sand.glm.lusi.summ)

# check the reference group for the factor, and change if you want
contrasts(env.qualit$LUSI.3000.impact)

# to change the factor's reference group: 
coef(manyglm(sand.mvabund.reduced ~ C(env.qualit$LUSI.3000.impact, base = 3), 
             family = "negative.binomial"))

## because in a glm it's hard to interpret the coefficients directly (glms use 
## link functions), better examine the predicted values
predict(sand.glm.lusi, type  =  "response")


anova.sand.glm.lusi <- anova.manyglm(sand.glm.lusi, 
                                     test = "LR",
                                     p.uni = "adjusted",
                                     blocks = tst.sites,
                                     nBoot = 50,
                                     show.time = "all")

print(anova.sand.glm.lusi$table)


# retrieve and plot species with the highest contribution to the patterns tested 
tst.lusi.uniSorted <- sort(anova.sand.glm.lusi$uni.test, index.return = TRUE, decreasing = TRUE)
plot(sand.mvabund.reduced ~ env.qualit$LUSI.3000.impact, var.subset = tst.lusi.uniSorted$ix[1:20]) 

# get the percentage of change in deviance due to these top ten/twenty => only about 
# 46%, so it's not right to focus on just them, we need more to explain the variability
sum(anova.sand.glm.lusi$uni.test[2, tst.lusi.uniSorted$ix[1:20]]) / sum(anova.sand.glm.lusi$uni.test[2, ])

# with 50 species, most of the change in deviance is explained (~78%)
sum(anova.sand.glm.lusi$uni.test[2, tst.lusi.uniSorted$ix[1:50]]) / sum(anova.sand.glm.lusi$uni.test[2, ])

############################################################################################################################
## try to recreate the plots, because originals using base graphics are kinda ugly
plot_mvabund <- function(abnd.df, gr.factor) {
  ## plots transformed abundance of species contributing the most to the 
  ## pattern studied by GLMs in mvabund (revision of plotting function from 
  ## package to make prettier plots).
  
  library(reshape2)
  library(ggplot2)
  
  # convert the df to long format -> easier for ggplot2 to handle  
  abnd.df.melted <- melt(abnd.df, id.vars = gr.factor) 
  
  # this doesn't keep row names, which here denote stations/replicates (does keep
  # row order, so easy to recover, if needed)
  
  ## use custom transformation for abundance values (log(y/min + 1)) - as in the 
  ## mvabund package
  # calculate the minimum non-0 value in the dataset 
  min.val <- min(abnd.df.melted$value[abnd.df.melted$value > 0]) 
  abnd.df.melted$value.tr <- log(abnd.df.melted$value / min.val + 1)
  
  ggplot(abnd.df.melted, aes_string(x = "variable", 
                                    y = "value.tr", 
                                    colour = gr.factor, 
                                    shape = gr.factor)) + 
    geom_point(size = 2) + 
    # reverse the order of x axis, so highest-contributing species are on top
    scale_x_discrete(name = "", limits = rev(levels(tst.df.melted$variable))) + 
    scale_y_continuous(name = "Abundance (log(y/min + 1))") + 
    coord_flip() + # put species on y axis - easier to read
    theme_bw()
}

## other kinds of raw data plots...


############################################################################################################################

# extract the data for plotting from the mvabund object - columns (species) now 
# sorted according to their contribution to explain the pattern tested -> subsetted, 
# because otherwise unreadable
tst.df <- as.data.frame(sand.mvabund.reduced[, tst.lusi.uniSorted$ix[1:50]])

# row order has not changed, so add the factor whose levels (= pattern) we tested 
tst.df$LUSI <- env.qualit$LUSI.3000.impact   

# plot & save (exploratory; will need subsetting to make it more presentable)
pdf(file.path(figs.dir, "explor_log-abnd_lusi_sand.pdf"),
    useDingbats = FALSE)
plot_mvabund(tst.df, "LUSI")
dev.off()

rm(tst.df)

## other simple abundance plots







## fit some other factor? - gravel/sand? sorting? O2 saturation? - COMBINATION OF FACTORS, MAYBE IDed FROM MDS ORDISURF? 


## see which species respond differently to different environmental parameters
## (= fit single predictive model for all species at all sites, but w/o attempting
## to explain the different responses using traits - the species ID is used in 
## place of a traits matrix). 
## Use a LASSO penalty (by setting the term method) when applying, to simplify 
## interpretation: does automatic model selection, setting to zero any interaction
## coefficients that don’t help reduce BIC.


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