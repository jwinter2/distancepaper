library(tidyverse)
library(arrow)
library(geosphere)
library(sf)

# Francois
setwd("~/Documents/Codes/distancepaper") # path to github repository

# Jordan
setwd("~/Downloads/distancepaper/")

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
      date < as.Date("2017-06-05") & cruise == "Gradients 2" ~ NaN, # contaminated salinity during first half the cruise
      c(abs(diff(salinity)),0) > 0.5 ~ NaN, # sporadic drops in salinity during TN271 cruise
      TRUE ~ salinity),
    temp = case_when(temp > 32 ~ NaN,
      TRUE ~ temp))

#-------------------
# Phytoplankton data
#-------------------
seaflow_psd <- arrow::read_parquet("data/PSD_TransitionZone_2022-05-24.parquet") %>%
  filter(pop != "croco") # remove this population as it is only found in low abundance in too few cruises

# Note for conversion from carbon per cell to equivalent spherical diameter
# Menden-Deuer, S. & Lessard, E. J. Carbon to volume relationships for dinoflagellates, diatoms, and other protist plankton. Limnol. Oceanogr. 45, 569–579 (2000).
d <- 0.261; e <- 0.860 # < 3000 µm3 

seaflow <- seaflow_psd %>% 
  group_by(cruise, date, pop) %>%
  dplyr::summarise(
    lat = unique(lat, na.rm = TRUE),
    lon = unique(lon, na.rm = TRUE),
    n_per_uL = sum(n_per_uL, na.rm = TRUE), # calculate cell abundance per population
    c_per_uL = sum(c_per_uL, na.rm = TRUE)) %>% # calculate carbon biomass per population
    mutate(qc = c_per_uL / n_per_uL, # calculate carbon per cell
      diam = round(2*(3/(4*base::pi)*(qc/d)^(1/e))^(1/3),5)) %>% # convert carbon per cell to equivalent spherical diameter *
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



#--------------------------------------
# A. Define the boundaries of the NPSG
#--------------------------------------

### Calculate distance along cruisetrack
# change longitude coordinate system for calculating distance
meta <- meta %>% 
  filter(!is.na(lon)) %>%
  mutate(lon = case_when(lon <= 0 ~ lon + 360,
    TRUE ~ lon)) %>%
  filter(lon > 100) # bad GPS coordinates


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
resolution <- 50 # km (i.e data binned every 10 km)
d <- range(meta_distance$distance)

meta_distance_binned <- meta_distance %>%
  group_by(
    cruise, 
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
  group_by(cruise) %>%
  dplyr::summarize(
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
    group_by(cruise) %>% 
    summarize(sal_front = max(sal_front, na.rm = TRUE),
      sal_front_up = max(sal_front_up, na.rm = TRUE),
      sal_front_down = max(sal_front_down, na.rm = TRUE)) 

 
### plot locations of all abrupt changes of salinity
# meta_distance_binned %>%  ggplot() + 
#   geom_point(aes(distance, salinity)) + 
#   geom_point(aes(distance, smooth_salinity), col = 2) + 
#   geom_point(aes(distance, smooth_salinity_down), col = 3) + 
#   geom_point(aes(distance, smooth_salinity_up), col = 4) + 
#   geom_hline(data = tz, aes( yintercept = sal_front), col = 2) + 
#   geom_hline(data = tz, aes( yintercept = sal_front_down), col = 3) + 
#   geom_hline(data = tz, aes( yintercept = sal_front_up), col = 4) + 
#   facet_wrap(. ~ cruise, scale = "free_x")


### Identify locations inside or outside the gyre (0 = gyre; 1 = outside gyre)
# Note: using the unbinned data this time
meta_gyre <- left_join(meta, tz) %>%
  group_by(cruise) %>%
  mutate(raw_gyre = case_when(salinity > sal_front ~ 0, 
      TRUE ~ 1),
    gyre = case_when(lon > 190 & lat < 24 & raw_gyre == 1 ~ 0,
      lon > 190 & lat < 11 & raw_gyre == 0 ~ 1,
      lon < 220 & lat < 33 & lat > 21 & raw_gyre == 1 ~ 0,
      lon < 190 & raw_gyre == 1 ~ 0,
      TRUE ~ raw_gyre)) %>%
  mutate(raw_gyre_down = case_when(salinity > sal_front_down ~ 0, 
      TRUE ~ 1),
    gyre_down = case_when(lon > 190 & lat < 24 & raw_gyre_down == 1 ~ 0,
      lon > 190 & lat < 11 & raw_gyre_down == 0 ~ 1,
      lon < 220 & lat < 33 & lat > 21 & raw_gyre_down == 1 ~ 0,
      lon < 190 & raw_gyre_down == 1 ~ 0,
      TRUE ~ raw_gyre_down)) %>%
  mutate(raw_gyre_up = case_when(salinity > sal_front_up ~ 0, 
      TRUE ~ 1),
    gyre_up = case_when(lon > 190 & lat < 24 & raw_gyre_up == 1 ~ 0,
      lon > 190 & lat < 11 & raw_gyre_up == 0 ~ 1,
      lon < 220 & lat < 33 & lat > 21 & raw_gyre_up == 1 ~ 0,
      lon < 190 & raw_gyre_up == 1 ~ 0,
      TRUE ~ raw_gyre_up)) %>%
  select(!c(raw_gyre, raw_gyre_up, raw_gyre_down, # remove unnecessary columns
    sal_front, sal_front_up, sal_front_down))


### calculate distance (km) from boundaries of the NPSG for each cruise
# location of boundaries
id <- which(diff(meta_gyre$gyre) != 0)
boundaries <- meta_gyre[id,] %>% 
  group_by(date) %>% 
  summarize(lat = mean(lat), 
    lon = mean(lon),
    cruise = unique(cruise))

id_down <- which(diff(meta_gyre$gyre_down) != 0)
boundaries_down <- meta_gyre[id_down,] %>% 
  group_by(date) %>% 
  summarize(lat = mean(lat), 
    lon = mean(lon),
    cruise = unique(cruise))

id_up <- which(diff(meta_gyre$gyre_up) != 0)
boundaries_up <- meta_gyre[id_up,] %>% 
  group_by(date) %>% 
  summarize(lat = mean(lat), 
    lon = mean(lon),
    cruise = unique(cruise))

### calculate distance from the front
sf_tz <- sf::st_as_sf(boundaries, coords = c("lon", "lat"), dim = 'XY', remove = FALSE, crs = 4326)
sf_meta <- sf::st_as_sf(meta_gyre, coords = c("lon", "lat"), dim = 'XY', remove = FALSE, crs = 4326)

meta_gyre_d <- tibble()
for(c in unique(sf_meta$cruise)) {
  sf_tz_g <- sf_tz %>% filter(cruise == c)
  sf_meta_g <- sf_meta %>% filter(cruise == c)
  meta_g <- meta_gyre %>% filter(cruise == c) %>%
    mutate(distance = apply(sf::st_distance(sf_meta_g, sf_tz), 1, function(x) min(x) / 1000)) %>%
    select(!c(gyre_down, gyre_up))
  meta_gyre_d <- bind_rows(meta_gyre_d, meta_g)
}

# calculate uncertainties in front location
sf_tz_down <- sf::st_as_sf(boundaries_down, coords = c("lon", "lat"), dim = 'XY', remove = FALSE, crs = 4326)
sf_tz_up <- sf::st_as_sf(boundaries_up, coords = c("lon", "lat"), dim = 'XY', remove = FALSE, crs = 4326)
boundaries <- boundaries %>% mutate(up = apply(sf::st_distance(sf_tz, sf_tz_up), 1, function(x) min(x) / 1000), 
  down = apply(sf::st_distance(sf_tz, sf_tz_down), 1, function(x) min(x) / 1000)) 

# transform to negative numbers when outside gyre
meta_gyre_d  <- meta_gyre_d %>% 
  mutate(distance = case_when(gyre == 1 ~ distance,
                              TRUE ~ - distance),
          gyre = case_when(gyre == 0 ~ "inside",
                           TRUE ~ "outside"))

### plot gyre
meta_gyre_d <- meta_gyre_d %>% filter(cruise != "SR1917" & cruise != "TN271")
g <- plot_geo(meta_gyre_d, lat = ~ lat, lon = ~ lon, color = ~ gyre, mode = "scatter", colors = c(viridis(4)[1],viridis(4)[3])) %>% layout(geo = geo)
#plotly_IMAGE(g, format = "png", out_file = "figures/cruise-track.png", width = 1000, height = 1000)















#---------------------------
# B. Calculate carbon growth
#---------------------------
### calculate rate of increase in carbon quotas during daylight (~ net carbon fixation)
meta_gyre_d <- meta_gyre_d %>% mutate(daytime = case_when(par > 10 ~ 1, TRUE ~ 0), # find daytime based on PAR values
  time_local = as.Date(date - 10 * 60 * 60)) %>%
  group_by(cruise, pop, daytime, time_local) %>%
  mutate(daily_growth = 60 * 60 * lm(qc ~ date)$coefficients[2] / mean(qc),
         growth_stderror = as.numeric(60 * 60 * broom::tidy(lm(qc ~ date))[2,3]),
         growth_pvalue = as.numeric(broom::tidy(lm(qc ~ date))[2,5]),
         period = as.numeric(diff(range(date)))) # period (in hours) of daylight

### curation of growth rate
meta_gyre_d <- meta_gyre_d %>%
  mutate(daily_growth = case_when( 
  daytime == 0 ~ NaN, # remove `growth` during nighttime
  growth_pvalue >= 0.01 ~ NaN, # NA's if p-value of grwoth rae is less than 0.01
  period < 6 ~ NaN, # NA's if daylight is less than 6 hours apart
  daily_growth < 0 ~ NaN, 
  TRUE ~ daily_growth)) %>%
  ungroup(time_local, daytime) %>%
  select(!c(time_local, daytime, period, growth_pvalue))

### plot daily growth estimates foreach cruise over time
meta_gyre_d %>% ggplot(aes(distance, daily_growth, col = pop)) +
  geom_pointrange(aes(ymin = daily_growth - growth_stderror, ymax = daily_growth + growth_stderror)) +
  facet_wrap(. ~ cruise, scale = "free_x") +
  theme_bw()

#--------------------------
# c. PLOTTING OVER DISTANCE
#--------------------------
### Calculate mean and sd over binned distance from the edges of the NPSG
res <- 100 # km (i.e data binned every 100 km)

d <- range(meta_gyre_d$distance)
data_figure <- meta_gyre_d %>%
  group_by(cruise, pop, 
    distance = cut(distance, seq(d[1],d[2], by = res), 
      labels = seq(d[1],d[2], by = res)[-1])) %>%
  dplyr::summarise_all(list(mean = function(x) mean(x, na.rm = TRUE), 
    sd = function(x) sd(x, na.rm = TRUE))) %>%
mutate(distance = as.numeric(as.character(distance))) %>%
ungroup(cruise)


### add East / South / North categories
data_figure <- data_figure %>%
  mutate(region = case_when(cruise == "SR1917" & lon < 220 | cruise == "Gradients 1" | cruise == "Gradients 2" | cruise == "Gradients 3" | cruise == "KM1712" | cruise == "KM1713" ~ "Northwest",
    cruise == "SR1917" & lon > 220 | cruise == "KM1502" | cruise == "TN271" | cruise == "TN398" | cruise == "Gradients 4" & lat > 22 ~ "Northeast",
    cruise == "Gradients 4" & lat < 22 | cruise == "KM1923" ~ "Southeast",
    TRUE ~ "Southwest"))

###### FOR NOW: don't use SR1917 and TN271
data_figure <- data_figure %>% filter(cruise != "SR1917" & cruise != "TN271")

### uncertainties in front location
binning <- res # uncertainties due to binning
front_uncertainties <- boundaries %>% 
filter(cruise != "SR1917" & cruise != "TN271") %>%
group_by(cruise) %>%
  summarize_all(mean) %>%
  mutate(up = case_when(up < binning ~ binning, 
      TRUE ~ up),
    down = case_when(down < binning ~ - binning, 
    TRUE ~ - down))

### choose a parameter
para <- "n_per_uL"; ylab <- "abundance (cells / μL)"; name <- "abundance"
# or
para <- "c_per_uL"; ylab <- "biomass (μgC / L)"; name <- "biomass"
# or
para <- "diam"; ylab <- "equivalent spherical diameter (μm)"; name <- "diameter"
# or
para <- "daily_growth"; ylab <- "daily growth rate"; name <- "growthrate"


### set colors
pop_cols <- c("prochloro" = viridis(3)[1], "picoeuk" = viridis(3)[2], "synecho" = viridis(3)[3])
gyre_cols <- c("inside" = viridis(5)[1], "outside" = viridis(5)[4])


### plotting parameter over distance per cruise
fig1 <- data_figure %>%
 select(region, cruise, pop, distance, contains(para)) %>%
 rename(mean = contains("mean"), sd = contains("sd")) %>%
  ggplot(aes(distance, mean,  col = pop)) + 
    geom_line(aes(group = pop), lwd = 1) +  
    geom_linerange(aes(ymax = mean + sd, ymin = mean - sd)) +
    geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
    # geom_vline(xintercept = 0, lty = 2) +
    scale_color_manual(values = pop_cols) +
    facet_wrap(. ~ cruise, scale = "free_x") +
    theme_bw() +
    labs(y = ylab, x = "distance (km)")

### plotting parameter over distance per region
fig2 <- data_figure %>%
 select(region, cruise, pop, distance, contains(para)) %>%
 rename(mean = contains("mean"), sd = contains("sd")) %>%
  ggplot(aes(distance, mean,  col = pop)) + 
    geom_line(aes(group = cruise), lwd = 1) +  
    geom_linerange(aes(ymax = mean + sd, ymin = mean - sd)) +
    geom_vline(xintercept = 0, lty = 2) +
    #scale_y_continuous(trans='log10') +
    scale_color_manual(values=pop_cols) +
    facet_grid(pop ~ region, scale = "free") +
    theme_bw() +
    labs(y = ylab, x = "distance (km)")

### plotting histogram of parameter inside vs outside the gyre
fig3 <- data_figure %>%
  select(region, cruise, pop, gyre, contains(para)) %>%
  rename(mean = contains("mean")) %>%
  ggplot(aes(x = mean, color = gyre, fill = gyre)) + 
  geom_histogram(aes(y = ..count../sum(..count..)), alpha = 0.4, bins=20, position = "identity") +
  facet_grid(pop ~ ., scale = "fixed") +
  scale_color_manual(values=gyre_cols) +
  scale_fill_manual(values=gyre_cols) +
  theme_bw() +
  #scale_x_continuous(trans='log10') +
  labs(x = ylab)

### plotting histogram of parameter inside vs outside the gyre per cruise
fig4 <- data_figure %>%
  select(region, cruise, pop, gyre, contains(para)) %>%
  rename(mean = contains("mean")) %>%
  ggplot(aes(x = mean, color = gyre, fill = gyre)) + 
  geom_histogram(aes(y = ..count../sum(..count..)), alpha = 0.4, bins=20, position = "identity") +  
  facet_grid(pop ~ cruise, scale = "fixed") +
  scale_color_manual(values=gyre_cols) +
  scale_fill_manual(values=gyre_cols) +
  theme_bw() +
  labs(x = ylab)

### plotting histogram of parameter inside vs outside the gyre per region
if (para == "c_per_uL"){ # get rid of outliers for easier visualization
  data_figure <- subset(data_figure, c_per_uL_mean < 60)
}
fig5 <- meta_gyre_d %>%
  select(region, cruise, pop, gyre, contains(para)) %>%
  rename(param = contains(para)) %>%
  ggplot(aes(x = param, color = gyre, fill = gyre)) + 
  geom_histogram(aes(y = ..count../sum(..count..)), alpha = 0.4, bins = 30, position = "identity") +  
  facet_grid(region ~ pop, scale = "free_x") +
  scale_color_manual(values=gyre_cols) +
  scale_fill_manual(values=gyre_cols) +
  theme_bw() +
  labs(x = ylab)


### save plot
png(paste0("figures/",name,"-distance-cruise.png"), width = 2500, height = 1600, res = 200)
print(fig1)
dev.off()

png(paste0("figures/",name,"-distance-region.png"), width = 2500, height = 1200, res = 200)
print(fig2)
dev.off()

png(paste0("figures/",name,"-gyre-hist.png"), width = 2500, height = 1200, res = 200)
print(fig3)
dev.off()

png(paste0("figures/",name,"-gyre-hist-cruise.png"), width = 2500, height = 800, res = 200)
print(fig4)
dev.off()

png(paste0("figures/",name,"-gyre-hist-region.png"), width = 2500, height = 1200, res = 200)
print(fig5)
dev.off()


### plotting nutrients
data_nutr <- data_figure[, c("region", "cruise", "pop", "distance", "NO3_NO2_mean", "PO4_mean")]
data_nutr <- data_nutr %>%
  rename(NO3_NO2 = NO3_NO2_mean, PO4 = PO4_mean) %>%
  gather(nutrient, concentration, 5:6)

fig_nutr <- data_nutr %>%
  ggplot(aes(distance, concentration, col = nutrient)) + 
  geom_point() +
  #geom_line(aes(group = nutrient), lwd = 1) +
  #geom_linerange(aes(ymax = mean + sd, ymin = mean - sd)) +
  scale_color_manual(values=c(viridis(4)[1], viridis(4)[3])) +
  geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
  theme_bw() +
  facet_wrap(. ~ cruise, scale = "free") +
  labs(y = "nutrient concentration (μg / L)", x = "distance (km)")

png(paste0("figures/",name,"-nutr-cruise.png"), width = 2500, height = 1200, res = 200)
print(fig_nutr)
dev.off()

# plotting individually

para <- "NO3_NO2"; ylab <- "nitrate and nitrite concentration (μg / L)"; name <- "nitrate"
para <- "PO4"; ylab <- "phosphate concentration (μg / L)"; name <- "phosphate"
#para <- "SiO4"; ylab <- "silicate concentration (μg / L)"; name <- "silicate"

### plotting parameter over distance per cruise

fig6 <- data_figure %>%
  select(region, cruise, pop, distance, contains(para)) %>%
  rename(mean = contains("mean"), sd = contains("sd")) %>%
  ggplot(aes(distance, mean)) + 
  geom_point() +
  geom_linerange(aes(ymax = mean + sd, ymin = mean - sd)) +
  geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
  theme_bw() +
  facet_wrap(. ~ cruise, scale = "free") +
  labs(y = ylab, x = "distance (km)")

### plotting histogram of parameter inside vs outside the gyre

fig7 <- data_figure %>%
  select(region, cruise, pop, gyre, contains(para)) %>%
  rename(mean = contains("mean")) %>%
  ggplot(aes(x = mean, color = gyre, fill = gyre)) + 
  geom_density(aes(y = ..scaled..), alpha = 0.4, bins = 50, position = "identity") +
  scale_color_manual(values=gyre_cols) +
  scale_fill_manual(values=gyre_cols) +
  theme_bw() +
  labs(x = ylab)

### save plot

png(paste0("figures/",name,"-gyre-cruise.png"), width = 2500, height = 1200, res = 200)
print(fig6)
dev.off()
  
png(paste0("figures/",name,"-gyre-hist.png"), width = 2500, height = 1200, res = 200)
print(fig7)
dev.off()


### plotting nutrients against phytoplankton diameter

### take mean diameter per day
meta_gyre_d$day <- substr(meta_gyre_d$date, 1, 10)
diam_data <- meta_gyre_d %>%
  group_by(cruise, pop, day, gyre) %>%
  dplyr::summarise_at(c("diam", "c_per_uL", "daily_growth", "SiO4", "NO3_NO2", "PO4"), list(mean), na.rm=T)


fig8 <- diam_data %>%
  select(cruise, pop, gyre, diam, contains(para)) %>%
  rename(mean = contains(para)) %>%
  ggplot(aes(mean, diam, col = gyre)) + 
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm", formula = y~x) +
  scale_color_manual(values=gyre_cols) +
  facet_wrap(pop ~ ., scale = "free") +
  labs(x = ylab, y = "daily diameter (μm)")


### plotting nutrients against phytoplankton biomass

fig9 <- diam_data %>%
  select(cruise, pop, gyre, c_per_uL, contains(para)) %>%
  rename(mean = contains(para)) %>%
  ggplot(aes(mean, c_per_uL, col = gyre)) + 
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm", formula = y~x) +
  scale_color_manual(values=gyre_cols) +
  facet_wrap(pop ~ ., scale = "free") +
  labs(x = ylab, y = "biomass (μgC / L)")

### plotting nutrients against growth rate

fig10 <- diam_data %>%
  select(cruise, pop, gyre, daily_growth, contains(para)) %>%
  rename(mean = contains(para)) %>%
  ggplot(aes(mean, daily_growth, col = gyre)) + 
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm", formula = y~x) +
  scale_color_manual(values=gyre_cols) +
  facet_wrap(pop ~ ., scale = "free") +
  labs(x = ylab, y = "daily growth rate")


### save plots

png(paste0("figures/",name,"-diam.png"), width = 2500, height = 1200, res = 200)
print(fig8)
dev.off()

png(paste0("figures/",name,"-biomass.png"), width = 2500, height = 1200, res = 200)
print(fig9)
dev.off()

png(paste0("figures/",name,"-growthrate.png"), width = 2500, height = 1200, res = 200)
print(fig10)
dev.off()
