import_zoo_data <- function(data.dir = NULL, zoo.data, station.names, repl) {
  ## imports, cleans and prepares zoobenthic raw data for community descriptive 
  ## statistics and further analyses. 
  ## 
  ## Arguments: data.dir - name of the directory where the dataset is
  ##              located (defaults to the working directory if not explicitly 
  ##              specified). 
  ##            zoo.data - name of csv file containing the species x stations 
  ##              data frame or matrix.
  ##            stations - vector of station names, in the order they should be
  ##              sorted in the final cleaned dataset. 
  ##            repl - number of replicates per station.
  ## Returns: data.frame in the form stations x species (format for community 
  ##  analyses used by vegan and similar functions and packages). There are no 
  ##  row names (rather, the row names are consecutive numbers). First 3 columns
  ##  are factors containing the station names, letters indicating the replicate 
  ##  number of the sample, and the year of sampling, respectively. The rest of 
  ##  the columns should be numeric, corresponding to the abundance of each species.
  ## Dependencies: gdata, plyr - for recoding the factors and sorting the output.
  
  library(plyr)
  library(gdata)
  
  # construct the file path. if there is a dedicated data subdirectory, look for the
  # dataset there; otherwise look in the working directory.
  if (is.null(data.dir)) {
    zoo.data.file <- zoo.data
  } else {
    zoo.data.file <- file.path(data.dir, zoo.data) 
  }
  
  # read the input raw data from the specified .csv file
  zoo.abnd.raw <- read.csv(file = zoo.data.file, header = TRUE, row.names = 1) 
  
  # transpose the data (stations x species matrix needed for vegan and vegan-dependent 
  # functions)
  zoo.abnd.transposed <- as.data.frame(t(zoo.abnd.raw))
  
  # get the horrible row names (labels that contain the station name, year and
  # replicate number) and extract all relevant info from them
  temp.row.names <- row.names(zoo.abnd.transposed)
  
  # get the station names: all of them end with a ., so truncate each string to the 
  # first .
  stations <- gsub(pattern = "\\..*", replacement = "", x = temp.row.names)
  
  # get the years: they are contained between two ., so split each string on the ., 
  # get only the second element (the year)
  years <- sapply(strsplit(temp.row.names, "\\."), "[", 2)
  
  # make the replicates vector: all the replicates for a single sampling of a 
  # station are assumed to be consecutive (i.e. you won't find replicate "a" of 
  # st1 coming after replicate "a" of st2, etc.). So, the replicates vector is 
  # repeated until it is the length of the stations vector. NB! will only work if 
  # there is always the same number of replicates/station. 
  replicates <- rep_len(letters[1:repl], length.out = length(stations))
  
  # get rid of the station names as row names
  row.names(zoo.abnd.transposed) <- c(1:nrow(zoo.abnd.transposed))  
  
  # construct the final data table containing everything we need 
  zoo.abnd.final <- as.data.frame(cbind(stations, 
                                        replicates, 
                                        years, 
                                        zoo.abnd.transposed))  
  
  # check for missing values. There shouldn't be any, so no values > 0 should 
  # appear here. If there are missing values, give an error message to that effect 
  # and print the table to see where the problem is.
  if (sum(apply(zoo.abnd.final, 2, function(x) sum(is.na(x)))) != 0) {
    print(apply(zoo.abnd.final, 2, function(x) sum(is.na(x))))
    stop("There are missing values in the data!")
    
  } else {
    
    # reorder the stations as desired
    zoo.abnd.final$stations <- reorder.factor(zoo.abnd.final$stations, 
                                              new.order = station.names)
    
    # sort the final data frame by station then year
    zoo.abnd.final <- arrange(zoo.abnd.final, stations, years)
    
    return(zoo.abnd.final)
  }
}