get_weather_data_pws <- function(weather.station, start.year, final.year)
{
  # Get long-term historical weather data from Weather Underground (www.wunderground.com) for the personal weather station (PWS) of your choice, in csv format. Can be modified for any station, knowing the start date (year) of the records, and the station code. 
  
  # set the station code  
  station.code <- weather.station
  
  # set the desired period for which you want the data.
  start.year <- start.year # The chosen start year/first year for which there is data for the station on the site; needs to be checked on the site beforehand. 
  final.year <- final.year
  years <- as.list(start.year : final.year)
  
  # prepare site web address. Since a maximum period of one year can be displayed, we will require the data from Jan 1 to Dec 31 and loop through the desired period year by year.
  site.prefix <- "https://www.wunderground.com/weatherstation/WXDailyHistory.asp?ID="
  site.middle1 <- "&day=1&month=1&year="
  site.middle2 <- "&dayend=31&monthend=12&yearend=" 
  site.suffix <- "&graphspan=custom&format=1"

  weather.temp.list <- lapply(years, function(x) {
                                 # construct the actual website address (only the year changes in the link).
                                 weather.link <- paste(site.prefix, station.code, site.middle1, x, site.middle2, x, site.suffix, sep = "")
    
                                 # retrieve the weather data as a data frame, and add it to the list of data
                                current.data <- read.csv(url(weather.link), sep = ",")
                                return(current.data) 
                                })
  
  # some of the column names differ between years, which will cause a problem for
  # combining the data frames later; fix by using a common vector of names for all
  # data frames (years)
  final.names <- names(weather.temp.list[[length(weather.temp.list)]])
  weather.temp.list <- lapply(weather.temp.list, function(x) {
                                                      x <- unname(x)
                                                      names(x) <- final.names
                                                      return(x)
                                                 })
  
  # merge the data frames from the list into a single data frame
  weather.data <- do.call("rbind", weather.temp.list)
  
  # convert the date column - the first column, to date type
  colnames(weather.data)[1] <- "Date"
  weather.data$Date <- as.Date(weather.data$Date)
  
  # every other row comes as containing all NAs; remove those - they are unnecessary
  weather.data <- weather.data[rowSums(is.na(weather.data)) != ncol(weather.data), ] # if the entire row consists of NAs, the sum should be equal to the number of columns in the data frame and so the row should be dropped.
  
  # insert a column containing the weather station code
  weather.data$station.code <- station.code
  
  # write the final data to file
  #write.csv(weather.data, file = "weather-data.csv")
  
  # final data
  return(weather.data)
}
