## mvabund - GLM fitting for abundance - environmental data
library(mvabund)

## checking model assumptions
# 1. mean-variance assumption => choice of family parameter. Can be checked by 
# plotting residuals vs fits: if little pattern - the chosen mean-variance 
# assumption is plausible.
# another way: direct plotting (variance ~ mean), for each species w/n each factor
# level.
# 2. assumed relationship b/n mean abundance and environmental variables - link 
# function and formula.
# If quantitative variables included in the model -> if trend in size of residuals
# at different fitted values (e.g. U-shape,..) = violation of the log-linearity 
# assumption.
#  

## testing hypotheses about the community-environment association
# -> anova on manyglm result

# advantage of multivariate analysis - greater power to detect patterns when 
# analysing all species simultaneously than when looking for a pattern specifically
# in each species.

