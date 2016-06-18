## summarize community abundance by larger tax.group

# 1. calculate prop of each species in each row (= sample) : tax.prop <- num.abund.df/rowSums(num.abund.df)
# 2. transpose -> species as rows
# 3. add new group identifier (new column) + other factors (station, year)
# 4. aggregate data by tax.group (= sum of proportions of all sp in that group)
# 5. maybe also aggregate by station and year (mean)?


# summarize abundance by year then station
summary.abnd.st.yr <- ddply(zoo.abnd.sand, .(years, stations), colwise(mean, .cols = is.numeric))

# calculate the proportion of each tax. group
head(summary.abnd.st.yr)
tst <- as.data.frame(t(summary.abnd.st.yr[sapply(summary.abnd.st.yr, is.numeric)]))
tst <- tst[rowSums(tst) > 0, ] 

tst$group <- current.zoo.taxa$group

tst.fin <- ddply(tst, .(group), transform, sum.n = length(group))

