# Tiny script for plotting weather data (temperatures) from Weather Underground
# weather.sub - data frame of weather data, column (variable) names unchanged from source (subsetted to only the study period). 

library(ggplot2)
library(scales)

# plot mean daily air temperature for the study period
ggplot(weather.sub, aes(x = Date, Mean.TemperatureC, colour = Mean.TemperatureC)) + 
  geom_line() + 
  scale_colour_gradient(low = "blue", high = "red", name = "Mean temperature (C)") + 
  scale_x_date(breaks = date_breaks("4 months"), labels = date_format("%b-%Y")) + # by default uses the current locale time setting for the labels; if you want the labels to match a specific locale, it can be set through Sys.setlocale(category = "LC_TIME", locale = your-choice-of-locale) 
  xlab("") + 
  ylab("Mean temperature (C)") +
  theme_bw()
