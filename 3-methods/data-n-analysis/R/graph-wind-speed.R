# Script to plot the mean daily wind speeds for the study period. Mostly identical to temperature script; consider making a small helper function to plot selected variables, only varying the colors and geoms.  
# weather.sub - data frame of weather data from Weather Underground (subsetted to study period only). 

library(ggplot2)
library(scales)

ggplot(weather.sub, aes(x = Date, Mean.Wind.SpeedKm.h, colour = Mean.Wind.SpeedKm.h)) + # color the points according to their value
  geom_point(alpha = 0.6, position = position_jitter(width = 3)) + # make the points transparent to be able to see them where they overlap. Won't work if saving to eps (doesn't support transparency) 
  scale_colour_gradient(name = "Mean wind speed (km/h)") + # set a nice name for the legend
  geom_smooth(method = "loess", size = 1) + # add a trendline, to make the midpoint of the data spread more easily visible
  scale_x_date(breaks = date_breaks("4 months"), labels = date_format("%b-%Y")) +  # by default uses the current locale time setting for the labels; if you want the labels to match a specific locale, it can be set through Sys.setlocale(category = "LC_TIME", locale = your-choice-of-locale) 
  xlab("") + 
  ylab("Mean wind speed (km/h)") +
  theme_bw() # make the background white (default is grey, gridded white)
