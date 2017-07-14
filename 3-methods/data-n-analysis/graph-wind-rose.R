##### Script to plot wind direction and speed as a wind rose in ggplot. 
##### Source & inspiration: Andy Clifton's (and the others) answers here:
##### http://stackoverflow.com/questions/17266780/wind-rose-with-ggplot-r

library(tidyverse)
library(forcats)

library(ggplot2)
library(viridis)

# plotting function
plot.windrose <- function(wind.data, wind.speed, wind.direction, 
                          speed.resolution = 2, 
                          direction.resolution = 22.5, 
                          speed.min = 2, 
                          speed.max = 50, 
                          speed.seq = NULL, 
                          palette = "YlGnBu", 
                          countmax = NA, 
                          debug = 0){
  
  # check the format of the data passed to the function
  if (is.numeric(wind.speed) & is.numeric(wind.direction)){
    # if the data is passed as vectors, construct a data frame 
    wind.data <- data.frame(wind.speed = wind.speed, wind.direction = wind.direction)
    speed <- "speed"
    direction <- "direction"
  }
  else if (exists("wind.data")){
    ## if the data is passed as a data frame, and the speed and direction are given as names. This is the preferred format.
  }
  
  ## tidy up the input data
  nb.in <- NROW(wind.data)
  dnu <- (is.na(wind.data[[wind.speed]]) | is.na(wind.data[[wind.direction]]))
  wind.data[[wind.speed]][dnu] <- NA
  wind.data[[wind.direction]][dnu] <- NA
  
  ## calculate the wind speed bins, if custom bins are not specified at function call
  if (missing(speed.seq)){
    speed.seq <- seq(speed.min, speed.max, speed.resolution)
  }
  else {
    if (debug > 0){
      cat("Using custom speed bins \n")
    }
  }
  
  ## get some information on number of bins, etc.
  nb.speed.seq <- length(speed.seq)
  
  if (max(wind.data[[wind.speed]], na.rm = TRUE) > speed.max){
    ## if the maximum wind speed in the dataset exceeds the default one (50), 
    ## use it instead, and adjust the bins accordingly.
    speed.breaks <- c(speed.seq, max(wind.data[[wind.speed]], na.rm = TRUE))
    speed.labels <- c(paste(c(speed.seq[1:nb.speed.seq - 1]), "-", c(speed.seq[2:nb.speed.seq])), 
                      paste(speed.max, "-", max(wind.data[[wind.speed]], na.rm = TRUE)))
  }
  else {
    speed.breaks <- speed.seq
    speed.labels <- paste(c(speed.seq[1:nb.speed.seq - 1]), "-", c(speed.seq[2:nb.speed.seq]))
  }
  
  wind.data$speed.binned <- cut(x = wind.data[[wind.speed]], 
                                breaks = speed.breaks, 
                                labels = speed.labels, 
                                ordered_result = TRUE)
  
  # figure out the direction bins
  direction.breaks <- c(0, 
                        seq(direction.resolution, 360 - direction.resolution, by = direction.resolution), 
                        360 + direction.resolution)
  direction.labels <- c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", 
                        "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")
  
  # assign each wind direction to a bin
  direction.binned <- cut(wind.data[[wind.direction]], 
                          breaks = direction.breaks, 
                          ordered_result = TRUE)
  levels(direction.binned) <- direction.labels
  wind.data$direction.binned <- direction.binned
  
  # run debug if required
  if (debug > 0){
    cat(direction.breaks, "\n")
    cat(direction.labels, "\n")
    cat(levels(direction.binned), "\n")
  }
  
  # create the plot
  windrose <- ggplot(data = wind.data, 
                     aes(x = direction.binned, fill = fct_rev(speed.binned))) +
                geom_bar() +
                scale_x_discrete(drop = FALSE, labels = waiver()) +
                coord_polar(start = -((direction.resolution)/360) * 2*pi) +
                theme_bw() +
                theme(axis.title.x = element_blank()) +  
                scale_fill_viridis(discrete = TRUE, direction = 1, name = "Wind speed (km/h)")
  
  # adjust axes if required
  if (!is.na(countmax)){
    windrose <- windrose + ylim(c(0, countmax))
  }
  
  # print the plot 
  print(windrose)
  
  # return the handle to the wind rose
  return(windrose)
}



#### plot the stupid wind rose, using 2012-2014 weather data
windrose.final.plot <- plot.windrose(wind.data = weather.sub,
                                     wind.speed = "Mean.Wind.SpeedKm.h",
                                     wind.direction = "WindDirDegrees",
                                     speed.seq = c(1, 6, 12, 20, 29, 39, 50)) + 
  # make separate subplots for each year (2012-2014)
  facet_wrap(~year) + 
  # reverse the order of the legend items
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme(legend.position = "bottom", 
        axis.text.x  = element_text(size = 6.25))

ggsave(file.path(figures.dir, "windrose_final.png"), 
       plot = windrose.final.plot, 
       dpi = 500, width = 16.5, units = "cm")
  
