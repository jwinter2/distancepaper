library(tidyverse)
library(arrow)
library(geosphere)
library(sf)
library(plotly)
library(viridis)
library(dplyr)
library(corrplot)

#-----------------
# For plotting map
#-----------------

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

setwd("/Users/ribalet/Documents/Codes/distancepaper")

#-------------------
# Environmental data
#-------------------
env <- read_csv("data/EnvironmentalData.csv") %>% arrange(date)
env$date <- as.POSIXct(env$date, format = "%m/%d/%y %H:%M", tz = "GMT")

### Clean-up data
# salinity and Temperature (due to freshwater contamination in ships' seawater supply)
env <- env %>%
  group_by(cruise) %>%
  mutate(
    salinity = case_when(salinity < 32 ~ NaN,
                         date < as.Date("2017-06-05") & cruise == "Gradients 2" ~ NaN, # contaminated salinity during first half the cruise
                         c(abs(diff(salinity)),0) > 0.5 ~ NaN, # sporadic drops in salinity during TN271 cruise
                         TRUE ~ salinity),
    temp = case_when(temp > 32 ~ NaN,
                     TRUE ~ temp))

# calculate daily-averaged par
env <- env %>%
  group_by(days = lubridate::floor_date(date, unit = "days")) %>%
  mutate(daily_par = mean(par)) %>%
  ungroup() %>%
  select(!days)

#-------------------
# Phytoplankton data
#-------------------
seaflow_psd <- arrow::read_parquet("data/PSD_TransitionZone_2022-05-24.parquet") %>% 
  mutate(n_per_uL = case_when(pop == "picoeuk" & cruise == "KM1923_751" ~ n_per_uL / 3, # adjust for virtual core correction
                              pop == "prochloro" & cruise == "KM1923_751" ~ n_per_uL * 2,
                              TRUE ~ n_per_uL),
         c_per_uL = case_when(pop == "picoeuk" & cruise == "KM1923_751" ~ c_per_uL / 3, # adjust for virtual core correction
                              pop == "prochloro" & cruise == "KM1923_751" ~ c_per_uL * 2 ,
                              TRUE ~ c_per_uL)) %>%
  filter(pop != "croco") %>% # remove this population as it is only found in low abundance in too few cruises
  group_by(pop) %>%
  mutate(pop = case_when(
    pop == "prochloro" ~ "Prochlorococcus",
    pop == "synecho" ~ "Synechococcus",
    pop == "picoeuk" & esd < 2 ~ "picoeukaryotes (< 2µm)",
    TRUE ~ "nanoeukaryotes (2-5µm)")) 


# Note for conversion from carbon per cell to equivalent spherical diameter
# Menden-Deuer, S. & Lessard, E. J. Carbon to volume relationships for dinoflagellates, diatoms, and other protist plankton. Limnol. Oceanogr. 45, 569–579 (2000).
d <- 0.261; e <- 0.860 # < 3000 µm3 

seaflow <- seaflow_psd %>% 
  group_by(cruise, date, pop) %>%
  dplyr::reframe(
    lat = unique(lat, na.rm = TRUE),
    lon = unique(lon, na.rm = TRUE),
    n_per_uL = sum(n_per_uL, na.rm = TRUE), # calculate cell abundance per population
    c_per_uL = sum(c_per_uL, na.rm = TRUE)) %>% # calculate carbon biomass per population
  mutate(qc = c_per_uL / n_per_uL, # calculate carbon per cell
         diam = round(2*(3/(4*base::pi)*(qc/d)^(1/e))^(1/3),5)) %>% # convert carbon per cell to equivalent spherical diameter *
  arrange(date) 

#-------------------
# Extract diel trend
#-------------------
decomp <-  function(data, day = 3){
  #print(unique(data[,"cruise"]))
  time <- unlist(list(zoo::rollapply(data[,"date"], width=24*day, function(x) unique(x), by.column = F)))
  diel <- unlist(list(zoo::rollapply(data[,"qc"], width=24*day, function(x) decompose(ts(x, frequency=24), type="mult")$seasonal, by.column= F))) -1
  quotas <- unlist(list(zoo::rollapply(data[,"qc"], width=24*day, function(x) decompose(ts(x, frequency=24), type="mult")$trend, by.column= F))) -1
  biomass <- unlist(list(zoo::rollapply(data[,"c_per_uL"], width=24*day, function(x) decompose(ts(x, frequency=24), type="mult")$trend, by.column= F))) -1
  abundance <- unlist(list(zoo::rollapply(data[,"n_per_uL"], width=24*day, function(x) decompose(ts(x, frequency=24), type="mult")$trend, by.column= F))) -1
  ts <- tibble(time, diel, quotas, biomass, abundance) %>%
    group_by(time) %>%
    summarize_all(function(x) mean(x, na.rm = T)) %>%
    ungroup()
  return(ts)
}

day <- 3
seaflow <-  seaflow %>%
  group_by(pop, cruise) %>%
  filter(length(unique(date)) > 24*day) %>%
  do(cbind(., decomp(., day))) %>%
  select(!time) %>%
  arrange(date) %>%
  ungroup(cruise) 


#----------------------------------------------------
# Combine environmental data with phytoplankton  data
#----------------------------------------------------
pre_meta <- full_join(seaflow %>% select(!cruise), env %>% select(!c(lat, lon)), by = "date") %>%
  group_by(cruise) %>%
  mutate(lat = zoo::na.approx(lat, na.rm = FALSE),
         lon = zoo::na.approx(lon, na.rm = FALSE))

# change longitude coordinate system 
pre_meta <- pre_meta  %>% 
  filter(!is.na(lon)) %>%
  mutate(lon = case_when(lon <= 0 ~ lon + 360,
                         TRUE ~ lon)) %>%
  filter(lon > 100) # bad GPS coordinates

# change Gradients cruise names to their official cruise names (R2R repository)
pre_meta  <- pre_meta  %>% mutate(cruise = case_when(cruise == "Gradients 1" ~ "KOK1606",
                                                     cruise == "Gradients 2" ~ "MGL1704",
                                                     cruise == "Gradients 3" ~ "KM1906",
                                                     cruise == "Gradients 4" & lat > 17 & lon > 215 ~ "TN397a",
                                                     cruise == "Gradients 4" & lat <= 17 & lon > 219.5 ~ "TN397b",
                                                     cruise == "Gradients 4" & lon <= 219.5 ~ "TN397c",
                                                     TRUE ~ cruise)) %>%
  filter(cruise != "SR1917") %>% # cruise outside of region of interest
  filter(cruise != "TN271") %>% # cruise outside of region of interest
  filter(cruise != "KM1502") # cruise with little variation


# Add cruise direction 
pre_meta  <- pre_meta  %>% mutate(direction = case_when(cruise == "KOK1606" | cruise == "MGL1704" | 
                                                        cruise == "KM1906" | cruise== "KM1712" |
                                                        cruise == "KM1713" ~ "North",
                                                        cruise == "TN398" | cruise == "TN397a"  ~ "East",
                                                        cruise == "TN397b" | cruise == "TN397c" | cruise == "KM1923" ~ "South"),
                                  direction = factor(direction, levels = c("North", "East", "South")),
                                  season = case_when(cruise == "KOK1606" | cruise == "MGL1704" | cruise == "KM1906" ~ "Spring",
                                                     cruise == "KM1712" | cruise == "KM1713" ~ "Summer",
                                                     cruise == "KM1923" | cruise == "TN397a" | cruise == "TN397c" | cruise == "TN397b" ~ "Fall",
                                                     cruise == "TN398" ~ "Winter"),
                                  season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter")),
                                  label = paste0(cruise, "/n(",direction, season,")"))

### plot cruise track
# plot_geo(pre_meta, lat = ~ lat, lon = ~ lon, color = ~ cruise, mode = "scatter") %>%  layout(geo = geo)



#-------------------------------
# Calculate cellular growth rate
#-------------------------------
### calculate rate of increase in carbon quotas during daylight (~ net carbon fixation)
meta <- pre_meta  %>% group_by(pop) %>%
  mutate(daytime = case_when(par > 10 ~ 1, 
                             TRUE ~ 0), # find daytime based on PAR values
         daynight = c(0,abs(diff(daytime))),
         day = cumsum(daynight)) %>%
  group_by(cruise, direction, pop, daytime, day) %>%
  mutate(n_obs = length(diel),
         daily_growth = case_when(daytime == 1 ~ diff(range(diel, na.rm= TRUE)),
                                  daytime == 0 ~ - diff(range(diel, na.rm= TRUE)))) %>%
  ungroup()


### curation of growth rate
meta <- meta %>%
  mutate(daily_growth = case_when( 
    daytime == 0 ~ NaN, # replace negative growth during nighttime by 0
    n_obs < 6 ~ NaN, # too few data point in a day
    TRUE ~ daily_growth)) %>%
  select(!c(day, daytime, daynight, n_obs))

### plot daily growth estimates for each cruise over time
# meta %>% ggplot(aes(date, daily_growth, col = pop)) +
#   geom_point(pch= 21) +
#   facet_wrap(. ~ cruise, scale = "free_x") +
#   theme_bw()






#--------------------------------------
# A. Define the boundaries of the NPSG
#--------------------------------------

### Calculate distance along cruisetrack
# ship's speed < 15 knots
# 1 km = 0.53996 knots (nautical mile / h) (top speed 14 knots ~ 26 km / h)
max_distance <- 14 / 0.53996 # km in 1 hour 

meta_distance <- tibble()
meta <- subset(meta, is.na(cruise) == F)

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
resolution <- 50 # km (i.e data binned every 10 km)
d <- range(meta_distance$distance)

meta_distance_binned <- meta_distance %>%
  group_by(
    cruise, direction, season, 
    distance = cut(distance, seq(d[1],d[2], by = resolution), 
      labels = seq(d[1],d[2], by = resolution)[-1])) %>%
  dplyr::summarise_all(function(x) mean(x, na.rm = TRUE)) %>%
  mutate(distance = as.numeric(as.character(distance)))


### plot salinity
# plot_geo(meta_distance, lat = ~ lat, lon = ~ lon, color = ~ salinity, mode = "scatter", colors = "Spectral") %>%  layout(geo = geo)


### Find abrupt changes in salinity (modified from Follett et al., 2021 Limnology & Oceanography 66 - 2442:2454, doi: 10.1002/lno.11763)
# calculate change of salinity over distance, using smooth salinity data
smooth_up <- 0.8 # upper smoothing parameter
smooth_down <- 0.6 # lower smoothing parameter
smooth <- mean(c(smooth_up,smooth_down)) # mean smoothing parameter

meta_distance_binned <- meta_distance_binned %>% 
  filter(!is.na(salinity) & !is.na(distance)) %>% # remove NA's for smoothing
  group_by(cruise, direction, season) %>%
  dplyr::reframe(
    date = date, lat = lat, lon = lon, distance = distance, salinity = salinity,
    smooth_salinity = smooth.spline(distance, salinity, all.knots= TRUE, spar = smooth)$y, # smooth data to remove noisy changes in salinity
    derivative_salinity = c(NA, abs(diff(smooth_salinity) / as.numeric(diff(distance)))), # calculate the rate of change of salinity over distance
    sal_front = case_when(derivative_salinity > quantile(derivative_salinity, 0.9, na.rm = TRUE) ~ smooth_salinity),  # find time of rapid changes in salinity
    
    smooth_salinity_up = smooth.spline(distance, salinity, all.knots= TRUE, spar = smooth_up)$y, # smooth data to remove noisy changes in salinity
    derivative_salinity_up = c(NA, abs(diff(smooth_salinity_up) / as.numeric(diff(distance)))), # calculate the rate of change of salinity over distance
    sal_front_up = case_when(derivative_salinity_up > quantile(derivative_salinity_up, 0.9, na.rm = TRUE) ~ smooth_salinity_up),   # find time of rapid changes in salinity
    
    smooth_salinity_down = smooth.spline(distance, salinity, all.knots= TRUE, spar = smooth_down)$y, # smooth data to remove noisy changes in salinity
    derivative_salinity_down = c(NA, abs(diff(smooth_salinity_down) / as.numeric(diff(distance)))), # calculate the rate of change of salinity over distance
    sal_front_down = case_when(derivative_salinity_down > quantile(derivative_salinity_down, 0.9, na.rm = TRUE) ~ smooth_salinity_down))   # find time of rapid changes in salinity
 
### Estimate salinity boundaries of NPSG's transition zone
# only look at the closest front to the NPSG's salinity (mean salinity ~ 35 psu) 
tz <- meta_distance_binned %>% 
    dplyr::group_by(cruise, direction, season) %>% 
    dplyr::reframe(sal_front = max(sal_front, na.rm = TRUE),
      sal_front_up = max(sal_front_up, na.rm = TRUE),
      sal_front_down = max(sal_front_down, na.rm = TRUE)) 

 
### plot locations of all abrupt changes of salinity
# meta_distance_binned %>%  ggplot() +
#   geom_point(aes(distance, salinity)) +
#   geom_point(aes(distance, smooth_salinity), col = "red") +
#   geom_point(aes(distance, smooth_salinity_down), col = "green") +
#   geom_point(aes(distance, smooth_salinity_up), col = "blue") +
#   geom_hline(data = tz, aes( yintercept = sal_front), col = "red") +
#   geom_hline(data = tz, aes( yintercept = sal_front_down), col = "green") +
#   geom_hline(data = tz, aes( yintercept = sal_front_up), col = "blue") +
#   facet_wrap(. ~ cruise, scale = "free_x")


### Identify locations inside or outside the gyre (0 = gyre; 1 = outside gyre)
# Note: using the unbinned data this time
meta_gyre <- left_join(meta, tz) %>%
  group_by(cruise, direction, season) %>%
  mutate(salinity = zoo::na.approx(salinity, na.rm = FALSE))

meta_gyre <- meta_gyre %>%
  mutate(raw_gyre = case_when(salinity > sal_front ~ 0,
                              TRUE ~ 1),
         gyre = case_when(
           lat < 6 & raw_gyre == 0 ~ 1,
           lat > 7 & lat < 23 & raw_gyre == 1 ~ 0,
           lon < 220 & lat < 33 & lat > 21 & raw_gyre == 1 ~ 0,
           TRUE ~ raw_gyre)) %>%
  select(!c(raw_gyre, sal_front, sal_front_down, sal_front_up)) # remove unnecessary columns


# plot gyre
# plot_geo(meta_gyre, lat = ~ lat, lon = ~ lon, color = ~ as.character(gyre), mode = "scatter") %>%  layout(geo = geo)


### calculate distance (km) from boundaries of the NPSG for each cruise
# location of boundaries
# the date is to make sure this is not finding a difference between different cruises
id <- which(diff(meta_gyre$gyre) != 0 & diff(meta_gyre$date) < 345600)
boundaries <- meta_gyre[id,] %>% 
  dplyr::group_by(date) %>% 
  dplyr::reframe(lat = mean(lat), 
    lon = mean(lon),
    cruise = unique(cruise),
    direction = unique(direction),
    season = unique(season),
    label = paste0(cruise, "\n(",direction, "/",season, ")"))

### calculate distance from the front
sf_tz <- sf::st_as_sf(boundaries, coords = c("lon", "lat"), dim = 'XY', remove = FALSE, crs = 4326)
sf_meta <- sf::st_as_sf(meta_gyre, coords = c("lon", "lat"), dim = 'XY', remove = FALSE, crs = 4326)

meta_gyre_d <- tibble()
for(c in unique(sf_meta$cruise)) {
  sf_tz_g <- sf_tz %>% filter(cruise == c)
  sf_meta_g <- sf_meta %>% filter(cruise == c)
  meta_g <- meta_gyre %>% filter(cruise == c) %>%
    mutate(distance = apply(sf::st_distance(sf_meta_g, sf_tz_g), 1, function(x) min(x) / 1000)) 
  meta_gyre_d <- bind_rows(meta_gyre_d, meta_g)
}


# transform to negative numbers when outside gyre
meta_gyre_d  <- meta_gyre_d %>% 
  mutate(distance = case_when(gyre == 1 ~ distance,
                              TRUE ~ - distance),
         gyre = case_when(
           gyre == 0 ~ "inside",
           gyre == 1 ~ "outside"),
         label = paste0(cruise, "\n(",direction, "/",season, ")"))

### remove data far outside the gyre
meta_gyre_d  <- meta_gyre_d %>% 
  filter(distance < 2000 & distance > -2000)



# plot gyre
#plot_geo(meta_gyre_d, lat = ~ lat, lon = ~ lon, color = ~ distance, mode = "scatter", colors =  colorRampPalette(c("blue2", "grey90", "red2"))(200)) %>%  layout(geo = geo)




#------------------------------
# b. BINNING data over DISTANCE
#------------------------------
meta_gyre_d <- read_csv("data/meta_gyre_d.csv")

# set colors
pop_cols <- c("Prochlorococcus" = rocket(7)[6], "Synechococcus" = rocket(7)[4], "picoeukaryotes (< 2µm)" = rocket(7)[3],  "nanoeukaryotes (2-5µm)" = rocket(7)[2])

### Calculate mean and sd over binned distance from the edges of the NPSG

binning <- 100 # kilometers 

data_figures <- meta_gyre_d  %>%
   group_by(cruise, pop, gyre, direction, season, label, distance = cut(distance, seq(-1000, 2000, by = binning), 
      labels = seq(-1000, 2000, by = binning)[-1]), .drop = T) %>%
  dplyr::summarise_all(list(mean = function(x) mean(x, na.rm = TRUE), 
    sd = function(x) sd(x, na.rm = TRUE))) %>%
  mutate(distance = as.numeric(as.character(distance)),
         pop = factor(pop, levels = c(names(pop_cols)[4], names(pop_cols)[3], names(pop_cols)[2], names(pop_cols)[1]))) %>%
  ungroup()

### uncertainties in front location due to binning
meta_gyre_d <- meta_gyre_d %>%
  mutate(gyre = case_when(
    abs(distance) < binning ~ "transition",
    TRUE ~ gyre))
data_figures <-data_figures %>%
  mutate(gyre = case_when(
    abs(distance) < binning ~ "transition",
    TRUE ~ gyre))

# plot cruise tracks
# plot_geo(meta_gyre_d, lat = ~ lat, lon = ~ lon, color = ~ gyre, mode = "scatter", colors = viridis(9, alpha = 1, begin = 0, end = 1, direction = 1)) %>%  layout(geo = geo)



# get proportion stats
prop <- data_figures %>% 
        group_by(distance, cruise, season) %>% 
        reframe(sum = sum(c_per_uL_mean))
prop <- merge(data_figures, prop)
prop <- prop %>% 
       group_by(pop, distance, cruise, season) %>% 
       reframe(prop_c = c_per_uL_mean / sum)




#----------------
# c. Main Figures
#----------------


### FIGURE 1

getPalette = colorRampPalette((RColorBrewer::brewer.pal(12, "Paired")))

# plot cruise track
fig1a <- meta_gyre_d %>%
  ggplot() +
  geom_path(aes(lon - 360, lat, color = cruise), lwd=2, show.legend = F) +
  coord_fixed(ratio = 1, xlim = c(-170, -110), ylim = c(-10, 60)) +
  borders("world", colour = "black", fill = "gray80") +
  theme_bw() +
  scale_color_manual(values = getPalette(10)) +
  theme(text = element_text(size = 20)) + 
  xlab("Longitude (ºW)") +
  ylab("Latitude (ºN)") +
  annotate("text", x=-140, y=50, label="KM1713") +
  annotate("text", x=-142, y=42, label="KM1712") +
  annotate("text", x=-166, y=45, label="KOK1606") +
  annotate("text", x=-166, y=40, label="MGL1704") +
  annotate("text", x=-166, y=35, label="KM1906") +
  annotate("text", x=-160, y=-1, label="KM1923") +
  annotate("text", x=-122, y=22, label="TN397a") +
  annotate("text", x=-133, y=10, label="TN397b") +
  annotate("text", x=-146, y=-2, label="TN397c") +
  annotate("text", x=-130, y=35, label="TN398")

# plot gyre boundaries
fig1b <- meta_gyre_d %>%
  ggplot() +
  geom_path(aes(lon - 360, lat, color = gyre, group = cruise), lwd = 2, show.legend = F) +
  coord_fixed(ratio = 1, xlim = c(-170, -110), ylim = c(-10, 60)) +
  borders("world", colour = "black", fill = "gray80") +
  theme_bw() +
  scale_color_manual(values = viridis(3)) +
  theme(text = element_text(size = 20)) + 
  xlab("Longitude (ºW)") +
  ylab("Latitude (ºN)")

# plot environmental variables

fig1c <- data_figures %>%
  filter(!is.na(salinity_mean)) %>%
  ggplot(aes(distance, salinity_mean, color = cruise)) +
  annotate("rect", xmin = -binning, xmax = binning, ymin = -Inf, ymax = Inf,  fill = "lightgrey") +
  geom_point(size = 2, show.legend = T) +
  geom_line(aes(col = cruise), lwd = 2) +
  geom_linerange(aes(ymin = salinity_mean - salinity_sd, ymax = salinity_mean + salinity_sd), lwd = 0.5) +
  theme_bw() +
  scale_color_manual(values = getPalette(10)) +
  theme(text = element_text(size = 20)) + 
  xlab("") +
  ylab("Salinity (psu)") 

fig1d <- data_figures %>%
  filter(!is.na(temp_mean)) %>%
  ggplot(aes(distance, temp_mean, color = cruise)) +
  annotate("rect", xmin = -binning, xmax = binning, ymin = -Inf, ymax = Inf,  fill = "lightgrey") +
  geom_point(size = 2,  show.legend = T) +
  geom_line(aes(col = cruise), lwd = 2) +
  geom_linerange(aes(ymin = temp_mean - temp_sd, ymax = temp_mean + temp_sd), lwd = 0.5) +
  theme_bw() +
  scale_color_manual(values = getPalette(10)) +
  theme(text = element_text(size = 20)) + 
  xlab("") +
  ylab("Temp (ºC)") 

fig1e <- data_figures %>%
  filter(!is.na(daily_par_mean)) %>%
  ggplot(aes(distance, daily_par_mean/100, color = cruise)) +
  annotate("rect", xmin = -binning, xmax = binning, ymin = -Inf, ymax = Inf,  fill = "lightgrey") +
  geom_point(size = 2,  show.legend = T) +
  geom_line(aes(col = cruise), lwd = 2) +
  geom_linerange(aes(ymin = daily_par_mean/100 - daily_par_sd/100, ymax = daily_par_mean/100 + daily_par_sd/100), lwd = 0.5) +
  theme_bw() +
  scale_color_manual(values = getPalette(10)) +
  theme(text = element_text(size = 20)) + 
  xlab("Distance (km)") +
  ylab(expression(paste("Light (µmolE m"^{-2},"s"^{-1},")")))
  

fig1f <- data_figures %>%
  filter(!is.na(NO3_NO2_mean)) %>%
  ggplot(aes(distance, NO3_NO2_mean, color = cruise)) +
  annotate("rect", xmin = -binning, xmax = binning, ymin = -Inf, ymax = Inf,  fill = "lightgrey") +
  geom_point(size = 2,  show.legend = T) +
  geom_line(aes(col = cruise), lwd = 2) +
  geom_linerange(aes(ymin = NO3_NO2_mean - NO3_NO2_sd, ymax = NO3_NO2_mean + NO3_NO2_sd), lwd = 0.5) +
  theme_bw() +
  scale_color_manual(values = getPalette(10)) +
  theme(text = element_text(size = 20)) + 
  xlab("") +
  ylab(expression(paste("DIN (µmol L"^{-1},")"))) +
  ylim(0,10)

fig1g <- data_figures %>%
  filter(!is.na(PO4_mean)) %>%
  ggplot(aes(distance, PO4_mean*10, color = cruise)) +
  annotate("rect", xmin = -binning, xmax = binning, ymin = -Inf, ymax = Inf,  fill = "lightgrey") +
  geom_point(size = 2,  show.legend = T) +
  geom_line(aes(col = cruise), lwd = 2) +
  geom_linerange(aes(ymin = PO4_mean*10 - PO4_sd*10, ymax = PO4_mean*10 + PO4_sd*10), lwd = 0.5) +
  theme_bw() +
  scale_color_manual(values = getPalette(10)) +
  theme(text = element_text(size = 20)) + 
  xlab("") +
  ylab(expression(paste("DIP (µmol L"^{-1},")"))) +
  ylim(0,10)

fig1h <- data_figures %>%
  filter(!is.na(MLD_mean)) %>%
  ggplot(aes(distance, MLD_mean, color = cruise)) +
  annotate("rect", xmin = -binning, xmax = binning, ymin = -Inf, ymax = Inf,  fill = "lightgrey") +
  geom_point(size = 2,  show.legend = T) +
  geom_line(aes(col = cruise), lwd = 2) +
  geom_linerange(aes(ymin = MLD_mean - MLD_sd, ymax = MLD_mean + MLD_sd), lwd = 0.5) +
  theme_bw() +
  scale_color_manual(values = getPalette(10)) +
  theme(text = element_text(size = 20)) + 
  xlab("Distance (km)") +
  ylab("Mixed Layer Depth (m)") 

png("figures/Figure_1.png", width = 2500, height = 3000, res = 200)
ggpubr::ggarrange(fig1a, fig1b, fig1c, fig1f, fig1d, fig1g,
                  ncol = 2, nrow = 3, 
                  common.legend = TRUE, legend = "bottom",
                  labels="AUTO") +
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm")) 
dev.off()


### FIGURE 2
# plot biomass over distance per cruise

# scale nutrient to biomass data
coeff <- 2

fig2 <- data_figures %>%
  ggplot(aes(distance, biomass_mean,   col = pop)) + 
  annotate("rect", xmin = -binning, xmax = binning, ymin = -Inf, ymax = Inf,  fill = "lightgrey") +
  geom_linerange(aes(ymin = biomass_mean - biomass_sd, ymax = biomass_mean + biomass_sd), lwd = 0.5,show.legend = FALSE) + 
  geom_line( aes(col = pop), lwd = 1, show.legend = TRUE) + 
  geom_point(aes(distance, NO3_NO2_mean * coeff), col = 1, pch = 16, size = 3, show.legend = FALSE) + 
  scale_color_manual(values = pop_cols, name = "population") +
  scale_y_continuous(name = expression(paste("Biomass (µgC L"^{-1},")")),
                     sec.axis = sec_axis(transform =~./coeff, 
                                         name = expression(paste("DIN (µmol L"^{-1},")")))) +
  facet_wrap(. ~ label, ncol = 5) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top") +
  xlab("Distance (km)")



png("figures/Figure_2.png", width = 2500, height = 1500, res = 200)
print(fig2)
dev.off()




### FIGURE 3
# Cellular Growth

fig3 <- data_figures %>%
  #filter(pop == "Prochlorococcus" | pop == "Synechococcus") %>%
  drop_na(daily_growth_mean) %>%
  ggplot(aes(distance, daily_growth_mean, group = pop)) +
  annotate("rect", xmin = -binning, xmax = binning, ymin = -Inf, ymax = Inf,  fill = "lightgrey") +
  geom_linerange(aes(ymin = daily_growth_mean - daily_growth_sd -  mean(daily_growth_sd, na.rm = TRUE), 
                     ymax = daily_growth_mean + daily_growth_sd +  mean(daily_growth_sd, na.rm = TRUE), col = pop), lwd= 0.5, show.legend = FALSE) +
  geom_line(aes(col= pop), lwd= 1, show.legend = TRUE) +
  scale_color_manual(values = pop_cols, name = "population") +
  facet_wrap(. ~ label, ncol = 5) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top") +
  labs(y = "Net cellular growth (1/d)", x = "Distance (km)")

png("figures/Figure_3.png", width = 2500, height = 1500, res = 200)
print(fig3)
dev.off()



### FIGURE 4
# correlation plot


# Function to process data for each direction
cols_to_keep <- c("light", "salinity", "temperature", "phosphate", "nitrate",
                  "growth Pro", "growth Syn",
                  "biomass Pro", "biomass Syn", "biomass nano", "biomass pico")
  
cor_all <-  meta_gyre_d  %>%
    select(pop, c_per_uL, diam, daily_growth, PO4, NO3_NO2, salinity, temp, daily_par) %>%
    pivot_wider(names_from = pop, values_from = c(c_per_uL, diam, daily_growth), values_fn = mean) %>%
    setNames(c("phosphate", "nitrate","salinity", "temperature", "light",
               "biomass Pro", "biomass Syn", "biomass nano", "biomass pico",
               "diam Pro", "diam Syn", "diam nano", "diam pico",
               "growth Pro", "growth Syn", "growth nano", "growth pico")) %>%
    select(all_of(cols_to_keep)) %>% {
      corr_data <- .
      cor_all <- cor(corr_data, use = "complete.obs")
      cor_all_p <- cor.mtest(corr_data, use = "complete.obs", conf.level = .99)
      
      # Adjust for multiple comparisons
      pAdj <- p.adjust(c(cor_all_p[[1]]), method = "BH")
      resAdj <- matrix(pAdj, ncol = dim(cor_all_p[[1]])[1])
      dimnames(resAdj) <- dimnames(cor_all_p$p)
      
      list(cor_matrix = cor_all, p_adj_matrix = resAdj)
    }





colors <- colorRampPalette(c("blue2", "grey90", "red2"))(200)

png("figures/Figure_4.png", width = 1200, height = 1200, res = 200)
par(mfrow=c(1,1))

corrplot(cor_all$cor_matrix, p.mat = cor_all$p_adj_matrix, 
         sig.level = 0.01, insig = "blank",
         diag = F,
         type = "lower", method = "color", addgrid.col = F, tl.col = "black", col = colors)



dev.off()
















# Function to process data for each direction
process_corr_data <- function(data, direction) {
  cols_to_keep <- c("light", "salinity", "temperature", "phosphate", "nitrate",
                    "growth Pro", "growth Syn",
                    "biomass Pro", "biomass Syn", "biomass nano", "biomass pico")
  
  meta_gyre_d  %>%
    filter(direction == {{ direction }}) %>%
    select(pop, c_per_uL, diam, daily_growth, PO4, NO3_NO2, salinity, temp, daily_par) %>%
    pivot_wider(names_from = pop, values_from = c(c_per_uL, diam, daily_growth), values_fn = mean) %>%
    setNames(c("phosphate", "nitrate","salinity", "temperature", "light",
               "biomass Pro", "biomass Syn", "biomass nano", "biomass pico",
               "diam Pro", "diam Syn", "diam nano", "diam pico",
               "growth Pro", "growth Syn", "growth nano", "growth pico")) %>%
    select(all_of(cols_to_keep)) %>% {
      corr_data <- .
      cor_all <- cor(corr_data, use = "complete.obs")
      cor_all_p <- cor.mtest(corr_data, use = "complete.obs", conf.level = .99)
      
      # Adjust for multiple comparisons
      pAdj <- p.adjust(c(cor_all_p[[1]]), method = "BH")
      resAdj <- matrix(pAdj, ncol = dim(cor_all_p[[1]])[1])
      dimnames(resAdj) <- dimnames(cor_all_p$p)
      
      list(cor_matrix = cor_all, p_adj_matrix = resAdj)
    }
}

# Process data for each direction
cor_all_north <- process_corr_data(data_figures, "North")
cor_all_east <- process_corr_data(data_figures, "East")
cor_all_south <- process_corr_data(data_figures, "South")



colors <- colorRampPalette(c("blue2", "grey90","red2"))(200)

png("figures/Figure_4.png", width = 2600, height = 1200, res = 200)

par(mfrow=c(1,3))
# North plot
corrplot(cor_all_north$cor_matrix, p.mat = cor_all_north$p_adj_matrix, 
         sig.level = 0.01, insig = "blank",
         diag = F,
         type = "lower", method = "color", addgrid.col = F, tl.col = "black", col = colors)
mtext("North", at=7, line = -38, cex=1)
mtext("a", at=-1, line = -2, cex=1, font=2)

# East plot
corrplot(cor_all_east$cor_matrix, p.mat = cor_all_east$p_adj_matrix, 
         sig.level = 0.01, insig = "blank",
         diag = F,
         type = "lower", method = "color", addgrid.col = F, tl.col = "black", col = colors)
mtext("East", at=7, line = -38, cex=1)
mtext("b", at=-1, line = -2, cex=1, font=2)

# South plot
corrplot(cor_all_south$cor_matrix, p.mat = cor_all_south$p_adj_matrix, 
         sig.level = 0.01, insig = "blank",
         diag = F,
         type = "lower", method = "color", addgrid.col = F, tl.col = "black", col = colors)
mtext("South", at=7, line = -38, cex=1)
mtext("c", at=-1, line = -2, cex=1, font=2)

dev.off()








#------------------------
# d. Supplemental Figures 
#------------------------

abund_pro <- data_figures %>%
     filter(distance > -1500 & pop == "Prochlorococcus") %>%
     ggplot(aes(distance, n_per_uL_mean,  col = pop, fill = pop)) + 
     geom_line(lwd = 1, position = "stack") + 
     geom_ribbon(aes(x=distance, y=n_per_uL_mean, ymin=n_per_uL_mean, ymax=n_per_uL_mean, group=pop, fill=pop), position="stack", alpha=0.5) +
     geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
     scale_fill_manual(values = pop_cols, name = "population") +
     scale_color_manual(values = pop_cols, guide = "none") +
     scale_y_continuous(name = "abundance (cells/μL)") +
     facet_wrap( direction ~ cruise, ncol = 5) +
     theme_bw(base_size = 13) +
     theme(legend.position = "top") +
     labs(y = "abundance (cells/μL)", x = "distance (km)")

abund_pico <- data_figures %>%
     filter(distance > -1500 & pop != "Prochlorococcus") %>%
     ggplot(aes(distance, n_per_uL_mean,  col = pop, fill = pop)) + 
     geom_line(lwd = 1, position = "stack") + 
     geom_ribbon(aes(x=distance, y=n_per_uL_mean, ymin=n_per_uL_mean, ymax=n_per_uL_mean, group=pop, fill=pop), position="stack", alpha=0.5) +
     geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
     scale_fill_manual(values = pop_cols, name = "population") +
     scale_color_manual(values = pop_cols, guide = "none") +
     scale_y_continuous(name = "abundance (cells/μL)") +
     facet_wrap( direction ~ cruise, ncol = 5) +
     theme_bw(base_size = 13) +
     theme(legend.position = "top") +
     labs(y = "abundance (cells/μL)", x = "distance (km)")

png("figures/Figure_abundance.png", width = 2500, height = 2500, res = 200)
ggpubr::ggarrange(abund_pro, abund_pico, ncol = 1, nrow = 2, 
                  labels="auto") +
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm")) 
dev.off()

cellsizepro <- data_figures %>%
  filter(distance > -1500) %>%
  filter(pop == "Prochlorococcus") %>%
  drop_na(diam_mean) %>%
  ggplot(aes(distance, diam_mean, group = pop)) +
  geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
  geom_line(aes(col= pop), lwd= 1, show.legend = TRUE) +
  geom_linerange(aes(ymin = diam_mean - diam_sd -  mean(diam_sd, na.rm = TRUE), 
                     ymax = diam_mean + diam_sd +  mean(diam_sd, na.rm = TRUE), col = pop), lwd= 0.5, show.legend = FALSE) +
  scale_color_manual(values = pop_cols, name = "population") +
  facet_wrap(direction ~ cruise, ncol = 5) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top") +
  labs(y = "Equivalent spherical diameter (μm)", x = "Distance (km)")


cellsizepico <- data_figures %>%
  filter(distance > -1500) %>%
  filter(pop != "Prochlorococcus") %>%
  drop_na(diam_mean) %>%
  ggplot(aes(distance, diam_mean, group = pop)) +
  geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
  geom_line(aes(col= pop), lwd= 1, show.legend = TRUE) +
  geom_linerange(aes(ymin = diam_mean - diam_sd -  mean(diam_sd, na.rm = TRUE), 
                     ymax = diam_mean + diam_sd +  mean(diam_sd, na.rm = TRUE), col = pop), lwd= 0.5, show.legend = FALSE) +
  scale_color_manual(values = pop_cols, name = "population") +
  facet_wrap(direction ~ cruise, ncol = 5) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top") +
  labs(y = "Equivalent spherical diameter (μm)", x = "Distance (km)")

png("figures/Figure_diameter.png", width = 2500, height = 2500, res = 200)
ggpubr::ggarrange(cellsizepro, cellsizepico, ncol = 1, nrow = 2, 
                  labels="auto") +
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm")) 
dev.off()

