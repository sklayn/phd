## LUSI

## REDO THIS CRAP WITHOUT A MILLION SEPARATE OBJECTS 

# import land use data & additional pressure scores
land.use.sand <- read.csv(file.path(data.dir, "lusi-test.csv"), header = TRUE)
other.pressures.sand <- read.csv(file.path(data.dir, "lusi-other-pressures.csv"), header = TRUE)
coast.corr <- read.csv(file.path(data.dir, "lusi-coast-corrections.csv"), header = TRUE)

# check for import errors
str(land.use.sand)
str(other.pressures.sand)
str(coast.corr)

# calculate LUSI
lusi.sand <- lusi(land.use.sand[-1], other.pressures.sand[-1], coast.corr[-1])

rm(land.use.sand, other.pressures.sand, coast.corr)
