library(tidyverse)
library(arrow)
library(geopshere)
library(sf)

setwd("~/Documents/Codes/distancepaper") # path to github repository

#-----------------
# For plotting map
#-----------------
library(plotly)
library(viridis)

geo <- list(
  showland = TRUE,
  showlakes = TRUE,
  showcountries = TRUE,
  showocean = TRUE,
  countrywidth = 0.5,
  landcolor = toRGB("grey90"),
  lakecolor = toRGB("white"),
  oceancolor = toRGB("white"),
  projection = list(
    type = 'orthographic',
    rotation = list(
      lon = -160,
      lat = 30,
      roll = 0
    )
  ),
  lonaxis = list(
    showgrid = TRUE,
    gridcolor = toRGB("gray40"),
    gridwidth = 0.5
  ),
  lataxis = list(
    showgrid = TRUE,
    gridcolor = toRGB("gray40"),
    gridwidth = 0.5
  )
)



#-------------------
# Environmental data
#-------------------
env <- read_csv("data/EnvironmentalData.csv") %>% arrange(date)

### Clean-up data
# salinity and Temperature (due to freshwater contamination in ships' seawater supply)
env <- env %>%
  group_by(cruise) %>%
  mutate(
    salinity = case_when(salinity < 32 ~ NaN,
      salinity < 33.5 & str_detect(cruise, pattern = "Gradients") ~ NaN,
      salinity < 34.5 & lat < 32 & cruise == "Gradients 2" ~ NaN,
      c(abs(diff(salinity)),0) > 0.5 ~ NaN, # sporadic drops in salinity in TN271
      TRUE ~ salinity),
    temp = case_when(temp > 32 ~ NaN,
      TRUE ~ temp))

#-------------------
# Phytoplankton data
#-------------------
seaflow_psd <- arrow::read_parquet("data/PSD_TransitionZone_2022-05-18.parquet") 

# Note for conversion from carbon per cell to equivalent spherical diameter
# Menden-Deuer, S. & Lessard, E. J. Carbon to volume relationships for dinoflagellates, diatoms, and other protist plankton. Limnol. Oceanogr. 45, 569–579 (2000).
d <- 0.261; e <- 0.860 # < 3000 µm3 

seaflow <- seaflow_psd %>% 
  group_by(cruise, date, pop) %>%
  dplyr::summarize(
    lat = unique(lat, na.rm = TRUE),
    lon = unique(lon, na.rm = TRUE),
    n_per_uL = sum(n_per_uL, na.rm = TRUE), # calculate cell abundance per population
    c_per_uL = sum(c_per_uL, na.rm = TRUE)) %>% # calculate carbon biomass per population
    mutate(qc = c_per_uL / n_per_uL, # calculate carbon per cell
      esd = round(2*(3/(4*base::pi)*(qc/d)^(1/e))^(1/3),5)) %>% # convert carbon per cell to equivalent spherical diameter *
    arrange(date) %>%
    ungroup(cruise) 

#----------------------------------------------------
# Combine environmental data with phytoplankton  data
#----------------------------------------------------
meta <- full_join(seaflow %>% select(!cruise), env %>% select(!c(lat, lon)), by = "date") %>%
        group_by(cruise) %>%
        mutate(lat = zoo::na.approx(lat, na.rm = FALSE),
              lon = zoo::na.approx(lon, na.rm = FALSE))


### plot cruise track
# plot_geo(meta, lat = ~ lat, lon = ~ lon, color = ~ cruise, mode = "scatter") %>%  layout(geo = geo)




### Calculate distance along cruisetrack
# change longitude coordinate system for calculating distance
meta <- meta %>% 
  filter(!is.na(lon)) %>%
  mutate(lon = case_when(lon < 0 ~ lon + 360,
                        TRUE ~ lon))

# ship's speed < 15 knots
# 1 km = 0.53996 knots (nautical mile / h) (top speed 14 knots ~ 26 km / h)
max_distance <- 14 / 0.53996 # km in 1 hour 

meta_distance <- tibble()
for(c in unique(meta$cruise)) {
  meta_c <- meta %>% filter(cruise == c)
  meta_c <- meta_c %>% 
    mutate(
      raw_distance = c(0, geosphere::distHaversine(as.matrix(meta_c[,c("lon","lat")]))/1000),
      capped_distance = case_when(raw_distance > max_distance ~ max_distance, # to prevent distance to exceed max distance
        TRUE ~ raw_distance), 
      distance = cumsum(capped_distance)) %>% 
    select(cruise, date, lat, lon, distance)

  meta_d <- left_join(meta_c, meta)
  meta_distance <- bind_rows(meta_distance, meta_d)
}



### calculate mean values over binned distance 
resolution <- 10 # km (i.e data binned every 10 km)
d <- range(meta_distance$distance)
n <- round(diff(d) / resolution)

meta_distance <- meta_distance %>%
  group_by(
    cruise, 
    distance = cut(distance, seq(d[1],d[2],length.out = n), 
      labels = seq(d[1],d[2],length.out = n -1))) %>%
  dplyr::summarize_all(function(x) mean(x, na.rm = TRUE)) %>%
  mutate(distance = as.numeric(as.character(distance)))








#----------------------------------
# Define the boundaries of the NPSG
#----------------------------------
# 
# Find latitude marked by abrupt changes in salinity (for reference, mean salinity for NPSG ~ 35 psu)
# (modified from Follett et al., 2021 Limnology & Oceanography 66 - 2442:2454, doi: 10.1002/lno.11763)

### plot salinity
# plot_geo(meta_distance, lat = ~ lat, lon = ~ lon, color = ~ salinity, mode = "scatter", colors = "Spectral") %>%  layout(geo = geo)

### calculate change iof salinity over distance, using smooth salinity data
smooth <- 0.7 # smoothing parameter
meta_distance <- meta_distance %>% 
  filter(!is.na(salinity) & !is.na(distance)) %>% # remove NA's for smoothing
  group_by(cruise) %>%
  dplyr::summarize(
    date = date, lat = lat, lon = lon, distance = distance, salinity = salinity,
    smooth_salinity = smooth.spline(distance, salinity, all.knots= TRUE, spar = smooth)$y, # smooth data to remove noisy changes in salinity
    derivative_salinity = c(NA, abs(diff(smooth_salinity) / as.numeric(diff(distance)))), # calculate the rate of change of salinity over distance
    tz = case_when(derivative_salinity > quantile(derivative_salinity, 0.9, na.rm = TRUE) ~ date))  # find time of rapid changes in salinity

# plot locations of all abrupt changes of salinity
# meta_distance %>%  ggplot() + 
#   geom_point(aes(distance, salinity)) + 
#   geom_point(aes(distance, smooth_salinity), col = 2) + 
#   geom_point(data = meta_distance %>% filter(!is.na(tz)), aes(distance, smooth_salinity), size = 3, col = 3) + 
#   facet_wrap(. ~ cruise, scale = "free_x")

### Estimate midpoint locations of the NPSG's limits for each cruise
tz <- meta_distance %>% 
  filter(!is.na(tz)) %>% 
  filter(lat > 24 | lat < 13) %>% # exclude salinity changes near Hawaiian Islands
  filter(lat < 50) %>% # exclude salinity fronts near Alaska coast
  filter(lon < 235) %>% # exclude salinity fronts near Pacific Northwest coast
  filter(lat > 3) %>% # exclude salinity fronts near Equatorial Pacifc  
  select(cruise, date, lat, lon, distance) %>% 
  group_by(cruise, distance = cut(distance, 10)) %>% 
  dplyr::summarize_all(function(x) median(x, na.rm = T)) %>% # select the midpoint across each salinity change
  filter(if(cruise[1] == "Gradients 4") TRUE else lat == max(lat)) %>%
  mutate(date = lubridate::floor_date(date, unit = "hour"),
    lat_tz = lat) %>%
  select(cruise, date, lat, lon, lat_tz) 
      

### plot location of NPSG boundaries
# plot_geo(meta_distance, lat = ~ lat, lon = ~ lon, color = ~ salinity, type = "scatter", mode = "scatter", colors = "Spectral") %>% 
# add_trace(data = tz, lat = ~ lat, lon = ~ lon, color = "salinity front", size  = 2) %>% layout(geo = geo)



#-----------------------------------------------
# Calculate distance from boundaries of the NPSG
#-----------------------------------------------

# Identify location inside or outside the gyre
meta_gyre <- left_join(meta, tz %>% select(cruise, date, lat_tz)) %>%
  group_by(cruise) %>%
  mutate(
    lat_tz = case_when(lat < 21 & lon < 215 ~ median(lat_tz, na.rm = TRUE),  # to deal with Gradients 4 cruise with multiple crossing of boundaries
      lat < 21 & lon > 215 ~ min(lat_tz, na.rm = TRUE),
      lat >= 21 ~ max(lat_tz, na.rm = TRUE)),
    gyre = case_when(cruise == "SR1917" & lon < 205 ~ "IN", 
      lat > 21 & lat > lat_tz ~ "OUT",
      lat < 21 & lat < lat_tz ~ "OUT",
      TRUE ~ "IN")) %>%
  select(!lat_tz)

# calculate distance (km) from boundaries of the NPSG for each cruise
sf_tz <- sf::st_as_sf(tz, coords = c("lon", "lat"), dim = 'XY', remove = FALSE, crs = 4326)
sf_meta <- sf::st_as_sf(meta_gyre, coords = c("lon", "lat"), dim = 'XY', remove = FALSE, crs = 4326)

meta_gyre_d <- tibble()
for(c in unique(sf_meta$cruise)) {
  sf_tz_g <- sf_tz %>% filter(cruise == c)
  sf_meta_g <- sf_meta %>% filter(cruise == c)
  meta_g <- meta_gyre %>% filter(cruise == c) %>%
    mutate(distance = apply(sf::st_distance(sf_meta_g, sf_tz), 1, function(x) min(x) / 1000)) 
  meta_gyre_d <- bind_rows(meta_gyre_d, meta_g)
}

# transform to negative numbers when outside gyre
meta_gyre_d  <- meta_gyre_d %>% 
  mutate(distance = case_when(gyre == "OUT" ~ distance,
                              TRUE ~ - distance))

### plot gyre
# plot_geo(meta_gyre_d, lat = ~ lat, lon = ~ lon, color = ~ gyre, mode = "scatter", colors = c("deeppink4","cyan4")) %>% layout(geo = geo)




#-----------------------
# PLOTTING OVER DISTANCE
#-----------------------

resolution <- 100 # km (i.e data binned every 100 km)

d <- range(meta_gyre_d$distance)
n <- round(diff(d) / resolution)

# by cruise
meta_gyre_d %>%
  group_by(
    cruise, 
    pop,
    distance = cut(distance, seq(d[1],d[2],length.out = n), 
      labels = seq(d[1],d[2],length.out = n -1))) %>%
  dplyr::summarize(
    value = mean(esd, na.rm = TRUE),
    sd = sd(esd, na.rm = TRUE)) %>%
  mutate(distance = as.numeric(as.character(distance))) %>%
  filter(pop != "croco") %>%
  ggplot(aes(distance, value, col = pop),) + 
    geom_line(lwd = 1) +  
    geom_pointrange(aes(ymax = value + sd, ymin = value - sd)) +
    geom_vline(xintercept = 0, lty = 2) +
    #scale_y_continuous(trans='log10') +
    facet_wrap(. ~  cruise, scale = "free") +
    theme_bw()


# by population
meta_gyre_d %>%
  group_by(
    cruise, 
    pop,
    distance = cut(distance, seq(d[1],d[2],length.out = n), 
      labels = seq(d[1],d[2],length.out = n -1))) %>%
  dplyr::summarize(
    value = mean(esd, na.rm = TRUE),
    sd = sd(esd, na.rm = TRUE)) %>%
  mutate(distance = as.numeric(as.character(distance))) %>%
  filter(pop != "croco") %>%
  ggplot(aes(distance, value, col = cruise),) + 
    geom_line(lwd = 1) +  
    geom_pointrange(aes(ymax = value + sd, ymin = value - sd)) +
    geom_vline(xintercept = 0, lty = 2) +
    #scale_y_continuous(trans='log10') +
    facet_wrap(. ~  pop, scale = "free") +
    theme_bw()
