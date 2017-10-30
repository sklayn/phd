import_zoo_data <- function(data.dir = NULL, zoo.data, meta.labels) {
  ## Imports, cleans and prepares zoobenthic raw data for community descriptive 
  ## statistics and further analyses. 
  ## 
  ## Arguments: data.dir - name of the directory where the dataset is
  ##              located (defaults to the working directory if not explicitly 
  ##              specified). 
  ##            zoo.data - name of csv file containing the species x stations 
  ##              data frame or matrix.
  ##            meta.labels - vector of labels for the metadata contained in the 
  ##             column names.
  ## NB Function depends on column names of the type 
  ##  "station.month.year.habitat.replicate" - contain the necessary metadata for 
  ##  each sampling.
  ## Output: data frame in the form stations x species (format for community 
  ##  analyses used by vegan and similar functions and packages). No row names! 
  ##  First 5 columns contain the metadata: station names, month and year of sampling, 
  ##  type of habitat, replicate code. The rest of the columns should be numeric, 
  ##  and should contain the numeric data on each species.
  ## Dependencies: tidyverse packages
  
  library(tidyverse)
  
  # construct the file path. if a dedicated data subdirectory is specified, 
  # look for the dataset there; otherwise look in the working directory.
  if (is.null(data.dir)) {
    zoo.data.file <- zoo.data
  } else {
    zoo.data.file <- file.path(data.dir, zoo.data) 
  }
  
  # read the input raw data from the specified .csv file
  zoo.abnd.raw <- read_csv(file = zoo.data.file) 
  
  # transpose the data (stations x species matrix needed for vegan and 
  # vegan-dependent functions). Done by gathering all columns except the first 
  # one (the species names) into a long data frame, then transposing back into 
  # a wide data frame.
  zoo.abnd.final <- zoo.abnd.raw %>% 
    gather(key = labels, value = values, 2:ncol(zoo.abnd.raw)) %>% 
    spread_(names(zoo.abnd.raw)[1], "values")

  # the labels that contain the station name, month, year, habitat type and 
  # replicate number are now in the variable (1st) column. Extract all relevant 
  # info from them. NB this leaves the numeric columns as character, but if needed,
  # can be converted later.
  zoo.abnd.final <- zoo.abnd.final %>% 
    separate(col = labels, into = meta.labels)

  # check for missing values. There shouldn't be any! If any are found, print the 
  # offending column names.
  if (sum(apply(zoo.abnd.final, 2, function(x) sum(is.na(x)))) != 0) {
    print(names(which(apply(zoo.abnd.final, 2, function(x) sum(is.na(x)) != 0))))
    stop("Missing values found!")
    
  } else {
    return(zoo.abnd.final)
  }
}
