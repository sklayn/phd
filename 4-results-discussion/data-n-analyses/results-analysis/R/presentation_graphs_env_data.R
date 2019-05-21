### This is a very quick script to make some illustrative graphs on the 
### environmental data at my study stations/sites. Doesn't deserve a separate 
### notebook, so leaving it here. 
### Going to use the pre-cleaned/calculated/compiled environmental data for each
### dataset. 

## import libraries
library(here)
library(tidyverse)
library(viridis)

## set ggplot theme to always black and white
theme_set(theme_bw())

## set working subdirectories
data.dir <- "data"
figures.dir <- "figs"
functions.dir <- "R"
save.dir <- "output"


##### Sozopol bay - 2012 #####
### import the sediment data 
(sediment.data.2012 <- read_csv(here(save.dir, "sediments_imputed_2012.csv")) %>% 
   ## convert station and habitat to factor, after recoding a little
   mutate(site = case_when(site == "Konski" ~ "K", 
                           site == "Ribka" ~ "R", 
                           site == "Gradina" ~ "G"), 
          habitat = case_when(habitat == "S" ~ "sand", 
                              habitat == "Z" ~ "Zostera")) %>% 
   mutate(site = factor(site, levels = c("K", "R", "G")), 
          habitat = factor(habitat, levels = c("sand", "Zostera")))  
)

## plot %TOM (selecting only that parameter and converting df to long)
(plot.tom.2012 <- ggplot(sediment.data.2012 %>% select(site:TOM) %>% gather(key = "var", value = "val", -c(site, habitat)), 
                         aes(x = site, y = val, fill = habitat)) + 
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(name = "", values = c("orange1", "light green")) + 
    labs(x = "Site", y = "%TOM") + 
    theme(legend.position = "top")
)

## save plot 
ggsave(here(figures.dir, "sediment_tom_2012.png"), 
       plot.tom.2012, 
       dpi = 300)

## plot sediment composition (grain sizes)
(plot.grain.size.2012 <- ggplot(sediment.data.2012 %>% select(site, habitat, gravel:silt_clay) %>% gather(key = "var", value = "val", -c(site, habitat)), 
                         aes(x = site, y = val, fill = var)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_viridis(name = "", discrete = T, option = "viridis", begin = 0.5) + 
    labs(x = "Site", y = "%") + 
    facet_wrap(~habitat) + 
    theme(legend.position = "top")
)

## save 
ggsave(here(figures.dir, "sediment_grain_size_2012.png"), 
       plot.grain.size.2012, 
       dpi = 300)


## plot mean grain sizes
(plot.mean.grain.size.2012 <- ggplot(sediment.data.2012 %>% select(site, habitat, mean_grain_size) %>% gather(key = "var", value = "val", -c(site, habitat)), 
                                aes(x = site, y = val, fill = habitat)) + 
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(name = "", values = c("orange", "light green")) + 
    labs(x = "Site", y = "Mean grain size (um)") + 
    theme(legend.position = "top")
)

## save 
ggsave(here(figures.dir, "sediment_mean_grain_size_2012.png"), 
       plot.mean.grain.size.2012, 
       dpi = 300)


##### Sand stations - Burgas Bay, 2013-2014 #####
### import the sediment data 
(sediment.data.smry.sand <- read_csv(here(save.dir, "sediments_imputed_sand.csv")) %>% 
   ## convert station to factor
   mutate(station = factor(station, levels = c("Kraimorie", "Chukalya", "Akin", "Sozopol", "Agalina", "Paraskeva"))) %>%
   ## average by station
   group_by(station) %>% 
   summarize_at(vars(TOM:Ni), list(~mean))
)


## plot %TOM (selecting only that parameter and converting df to long)
(plot.tom.sand <- ggplot(sediment.data.smry.sand %>% select(station, TOM), 
                         aes(x = station, y = TOM)) + 
    geom_bar(stat = "identity", fill = "orange") +
    labs(x = "Station", y = "%TOM") + 
    scale_x_discrete(labels = paste0("S", as.numeric(unique(sediment.data.smry.sand$station))))
)

## save plot 
ggsave(here(figures.dir, "sediment_tom_sand.png"), 
       plot.tom.sand, 
       dpi = 300)

## plot sediment composition (grain sizes)
(plot.grain.size.sand <- ggplot(sediment.data.smry.sand %>% select(station, gravel:silt_clay) %>% gather(key = "var", value = "val", -station), 
                                aes(x = station, y = val, fill = var)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_viridis(name = "", discrete = T, option = "viridis", begin = 0.5) + 
    labs(x = "Station", y = "%") +
    scale_x_discrete(labels = paste0("S", as.numeric(unique(sediment.data.smry.sand$station)))) +
    theme(legend.position = "top")
)

## save 
ggsave(here(figures.dir, "sediment_grain_size_sand.png"), 
       plot.grain.size.sand, 
       dpi = 300)


## plot mean grain sizes
(plot.mean.grain.size.sand <- ggplot(sediment.data.smry.sand %>% select(station, mean_grain_size), 
                                     aes(x = station, y = mean_grain_size)) + 
    geom_bar(stat = "identity", fill = "sky blue") +
    labs(x = "Station", y = "Mean grain size (um)") + 
    scale_x_discrete(labels = paste0("S", as.numeric(unique(sediment.data.smry.sand$station))))   
)

## save 
ggsave(here(figures.dir, "sediment_mean_grain_size_sand.png"), 
       plot.mean.grain.size.sand, 
       dpi = 300)



##### Seagrass stations - Burgas Bay, 2013-2014 #####
### import the sediment data 
(sediment.data.smry.zostera <- read_csv(here(save.dir, "sediments_imputed_seagrass.csv")) %>% 
   ## convert station to factor
   mutate(station = factor(station, levels = c("Poda", "Otmanli", "Vromos", "Gradina", "Ropotamo"))) %>%
   ## average by station
   group_by(station) %>% 
   summarize_at(vars(TOM:Ni), list(~mean))
)


## plot %TOM (selecting only that parameter and converting df to long)
(plot.tom.zostera <- ggplot(sediment.data.smry.zostera %>% select(station, TOM), 
                         aes(x = station, y = TOM)) + 
    geom_bar(stat = "identity", fill = "orange") +
    labs(x = "Station", y = "%TOM") + 
    scale_x_discrete(labels = paste0("Z", as.numeric(unique(sediment.data.smry.zostera$station))))
)

## save plot 
ggsave(here(figures.dir, "sediment_tom_zostera.png"), 
       plot.tom.zostera, 
       dpi = 300)

## plot sediment composition (grain sizes)
(plot.grain.size.zostera <- ggplot(sediment.data.smry.zostera %>% select(station, gravel:silt_clay) %>% gather(key = "var", value = "val", -station), 
                                aes(x = station, y = val, fill = var)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_viridis(name = "", discrete = T, option = "viridis", begin = 0.5) + 
    labs(x = "Station", y = "%") +
    scale_x_discrete(labels = paste0("Z", as.numeric(unique(sediment.data.smry.zostera$station)))) +
    theme(legend.position = "top")
)

## save 
ggsave(here(figures.dir, "sediment_grain_size_zostera.png"), 
       plot.grain.size.zostera, 
       dpi = 300)


## plot mean grain sizes
(plot.mean.grain.size.zostera <- ggplot(sediment.data.smry.zostera %>% select(station, mean_grain_size), 
                                     aes(x = station, y = mean_grain_size)) + 
    geom_bar(stat = "identity", fill = "sky blue") +
    labs(x = "Station", y = "Mean grain size (um)") + 
    scale_x_discrete(labels = paste0("Z", as.numeric(unique(sediment.data.smry.zostera$station))))   
)

## save 
ggsave(here(figures.dir, "sediment_mean_grain_size_zostera.png"), 
       plot.mean.grain.size.zostera, 
       dpi = 300)