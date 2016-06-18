get_weather_data_airport <- function(weather.station, start.year, final.year)
{
  # Get long-term historical weather data from Weather Underground (www.wunderground.com) for the airport weather station of your choice, in csv format. Can be modified for any airport station, knowing the start date (year) of the records, and the station code. 
  
  # set the station code  
  station.code <- weather.station
  
  # set the desired period for which you want the data.
  start.year <- start.year # The chosen start year/first year for which there is data for the station on the site; needs to be checked on the site beforehand. 
  final.year <- final.year
  years <- as.list(start.year : final.year)
  
  # prepare site web address. Since a maximum period of one year can be displayed, we will require the data from Jan 1 to Dec 31 and loop through the desired period year by year.
  site.prefix <- "https://www.wunderground.com/history/airport/"
  site.middle <- "/1/1/CustomHistory.html?dayend=31&monthend=12&yearend="
  site.suffix <- "&req_city=&req_state=&req_statename=&reqdb.zip=&reqdb.magic=&reqdb.wmo=&format=1"
  
  weather.temp.list <- lapply(years, function(x) {
                                 # construct the actual website address (only the year changes in the link).
                                 weather.link <- paste(site.prefix, station.code, "/", x, site.middle, x, site.suffix, sep = "")
    
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
  
  # the break character has come as part of the last column and caused it to be read in as a factor; we'll remove it and convert the column back to numeric type 
  colnames(weather.data)[23] <- "WindDirDegrees"
  weather.data$WindDirDegrees <- lapply(weather.data$WindDirDegrees, function(y) gsub("<br />", "", y)) 
  weather.data$WindDirDegrees <- as.numeric(weather.data$WindDirDegrees)
  
  # final data
  return(weather.data)
}
