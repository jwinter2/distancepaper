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
env<- env %>%
  group_by(days = lubridate::floor_date(date, unit = "days")) %>%
  mutate(daily_par = mean(par)) %>%
  ungroup() %>%
  select(!days)

#-------------------
# Phytoplankton data
#-------------------
seaflow_psd <- arrow::read_parquet("data/PSD_TransitionZone_2022-05-24.parquet") %>%
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
  dplyr::summarise(
    lat = unique(lat, na.rm = TRUE),
    lon = unique(lon, na.rm = TRUE),
    n_per_uL = sum(n_per_uL, na.rm = TRUE), # calculate cell abundance per population
    c_per_uL = sum(c_per_uL, na.rm = TRUE)) %>% # calculate carbon biomass per population
    mutate(qc = c_per_uL / n_per_uL, # calculate carbon per cell
      diam = round(2*(3/(4*base::pi)*(qc/d)^(1/e))^(1/3),5)) %>% # convert carbon per cell to equivalent spherical diameter *
    arrange(date) 

#------------------------------------
# Extract diel trend of carbon quotas
#------------------------------------
decomp <-  function(data, day = 3){
  #print(unique(data[,"cruise"]))
  time <- unlist(list(zoo::rollapply(data[,"date"], width=24*day, function(x) unique(x), by.column = F)))
  diel <- unlist(list(zoo::rollapply(data[,"qc"], width=24*day, function(x) decompose(ts(x, frequency=24), type="mult")$seasonal, by.column= F))) -1
  ts <- tibble(time, diel) %>%
    group_by(time) %>%
    summarize_all(list(mean = mean, sd = sd)) %>%
    ungroup()
  return(ts)
}

day <- 3
seaflow <-  seaflow %>%
  group_by(pop, cruise) %>%
  filter(length(unique(date)) > 24*day) %>%
  do(cbind(., decomp(., day))) %>%
  rename(qc_diel = mean,
         qc_diel_sd = sd) %>%
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
                                                cruise == "KM1713" ~ "north",
                                              cruise == "TN398" | cruise == "TN397a" |
                                                cruise == "TN397b" ~ "east",
                                              cruise == "TN397c" | cruise == "KM1923" ~ "south"),
           direction = factor(direction, levels = c("north", "east", "south")))

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
  mutate(n_obs = length(qc_diel),
         daily_growth = case_when(daytime == 1 ~ diff(range(qc_diel, na.rm= TRUE)),
                                  daytime == 0 ~ - diff(range(qc_diel, na.rm= TRUE))))
  

### curation of growth rate
meta <- meta %>%
  mutate(daily_growth = case_when( 
    daytime == 0 ~ NaN, # replace negative growth during nighttime by 0
    n_obs < 6 ~ NaN, # too few data point in a day
    TRUE ~ daily_growth)) %>%
  ungroup(daytime, day, cruise, direction, pop) %>%
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
    cruise, direction, 
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
  group_by(cruise, direction) %>%
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
    dplyr::group_by(cruise, direction) %>% 
    dplyr::summarize(sal_front = max(sal_front, na.rm = TRUE),
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
  group_by(cruise, direction) %>%
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
# the date is to make sure this is not finding a difference between different cruises
id <- which(diff(meta_gyre$gyre) != 0 & diff(meta_gyre$date) < 345600)
boundaries <- meta_gyre[id,] %>% 
  dplyr::group_by(date) %>% 
  dplyr::summarize(lat = mean(lat), 
    lon = mean(lon),
    cruise = unique(cruise),
    direction = unique(direction))

id_down <- which(diff(meta_gyre$gyre_down) != 0 & diff(meta_gyre$date) < 345600)
boundaries_down <- meta_gyre[id_down,] %>% 
  dplyr::group_by(date) %>% 
  dplyr::summarize(lat = mean(lat), 
    lon = mean(lon),
    cruise = unique(cruise),
    direction = unique(direction))

id_up <- which(diff(meta_gyre$gyre_up) != 0 & diff(meta_gyre$date) < 345600)
boundaries_up <- meta_gyre[id_up,] %>% 
  dplyr::group_by(date) %>% 
  dplyr::summarize(lat = mean(lat), 
    lon = mean(lon),
    cruise = unique(cruise),
    direction = unique(direction))


### calculate distance from the front
sf_tz <- sf::st_as_sf(boundaries, coords = c("lon", "lat"), dim = 'XY', remove = FALSE, crs = 4326)
sf_meta <- sf::st_as_sf(meta_gyre, coords = c("lon", "lat"), dim = 'XY', remove = FALSE, crs = 4326)

meta_gyre_d <- tibble()
for(c in unique(sf_meta$cruise)) {
  sf_tz_g <- sf_tz %>% filter(cruise == c)
  sf_meta_g <- sf_meta %>% filter(cruise == c)
  meta_g <- meta_gyre %>% filter(cruise == c) %>%
    mutate(distance = apply(sf::st_distance(sf_meta_g, sf_tz_g), 1, function(x) min(x) / 1000)) %>%
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

### remove data far outside the gyre
meta_gyre_d  <- meta_gyre_d %>% filter(distance < 2000)




#--------------------------
# b. BINNING OVER DISTANCE
#--------------------------
# set colors
pop_cols <- c("Prochlorococcus" = rocket(7)[6], "Synechococcus" = rocket(7)[4], "picoeukaryotes (< 2µm)" = rocket(7)[3],  "nanoeukaryotes (2-5µm)" = rocket(7)[2])
### Calculate mean and sd over binned distance from the edges of the NPSG

res <- 250 # km (i.e data binned every 100 km)

d <- range(meta_gyre_d$distance)

data_figures <- meta_gyre_d %>%
   group_by(cruise, pop, direction,
    distance = cut(distance, seq(d[1],d[2], by = res), 
      labels = seq(d[1],d[2], by = res)[-1])) %>%
  dplyr::summarise_all(list(mean = function(x) mean(x, na.rm = TRUE), 
    sd = function(x) sd(x, na.rm = TRUE))) %>%
  mutate(distance = as.numeric(as.character(distance)),
         pop = factor(pop, levels = c(names(pop_cols)[4], names(pop_cols)[3], names(pop_cols)[2], names(pop_cols)[1]))) %>%
  ungroup(cruise, direction)


### uncertainties in front location
binning <- res # uncertainties due to binning
front_uncertainties <- boundaries %>% 
  filter(cruise != "SR1917" & cruise != "TN271") %>%
  group_by(cruise, direction) %>%
  summarize_all(mean) %>%
  mutate(up = case_when(up < binning ~ binning, 
      TRUE ~ up),
    down = case_when(down < binning ~ - binning, 
    TRUE ~ - down))

for (cruise_name in unique(front_uncertainties$cruise)){
  small_unc <- front_uncertainties %>% filter(cruise == cruise_name)
  up <- which(meta_gyre_d$cruise == cruise_name & meta_gyre_d$distance <= small_unc$up
              & meta_gyre_d$distance >= 0)
  meta_gyre_d$gyre[up] <- "transition"
  down <- which(meta_gyre_d$cruise == cruise_name & meta_gyre_d$distance >= small_unc$down
              & meta_gyre_d$distance <= 0)
  meta_gyre_d$gyre[down] <- "transition"
}



# plot cruise tracks
g <- plot_geo(meta_gyre_d, lat = ~ lat, lon = ~ lon, color = ~ gyre, mode = "scatter", colors = viridis(9, alpha = 1, begin = 0, end = 1, direction = 1)) %>%  layout(geo = geo)
#g

# get proportion stats
prop <- data_figures %>% 
        group_by(distance, cruise) %>% 
        reframe(sum = sum(c_per_uL_mean))
prop <- merge(data_figures, prop)
prop <- prop %>% 
       group_by(pop, distance, cruise) %>% 
       reframe(prop_c = c_per_uL_mean / sum)

#----------------
# c. Main Figures
#----------------


### FIGURE 1

getPalette = colorRampPalette((RColorBrewer::brewer.pal(12, "Paired")))

# plot cruise track
fig1a <- meta_gyre_d %>%
  filter(distance > -2000) %>%
  ggplot() +
  geom_point(aes(lon - 360, lat, color = cruise), size=2, alpha = 0.7, show.legend = F) +
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
  filter(distance > -2000) %>%
  ggplot() +
  geom_point(aes(lon - 360, lat, color = gyre), size=2, alpha = 0.7, show.legend = F) +
  coord_fixed(ratio = 1, xlim = c(-170, -110), ylim = c(-10, 60)) +
  borders("world", colour = "black", fill = "gray80") +
  theme_bw() +
  scale_color_manual(values = viridis(3)) +
  theme(text = element_text(size = 20)) + 
  xlab("Longitude (ºW)") +
  ylab("Latitude (ºN)")

# plot environmental variables
fig1c <- data_figures %>%
  filter(distance > -1500) %>%
  filter(!is.na(NO3_NO2_mean)) %>%
  ggplot(aes(distance, NO3_NO2_mean, color = cruise)) +
  geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.1, inherit.aes = FALSE) +
  geom_point(size = 2,  show.legend = T) +
  geom_line(aes(col = cruise), na.rm = TRUE, lwd = 2) +
  geom_linerange(aes(ymin = NO3_NO2_mean - NO3_NO2_sd, ymax = NO3_NO2_mean + NO3_NO2_sd), lwd = 0.5) +
  theme_bw() +
  scale_color_manual(values = getPalette(10)) +
  theme(text = element_text(size = 20)) + 
  xlab("Distance (km)") +
  ylab("DIN (µmol/L)") 

fig1d <- data_figures %>%
  filter(distance > -1500) %>%
  ggplot(aes(distance, temp_mean, color = cruise)) +
  geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.1, inherit.aes = FALSE) +
  geom_point(size = 2,  show.legend = T) +
  geom_line(aes(col = cruise), lwd = 2) +
  geom_linerange(aes(ymin = temp_mean - temp_sd, ymax = temp_mean + temp_sd), lwd = 0.5) +
  theme_bw() +
  scale_color_manual(values = getPalette(10)) +
  theme(text = element_text(size = 20)) + 
  xlab("Distance (km)") +
  ylab("Temperature (ºC)") 


png("figures/Figure_1.png", width=12, height=12, unit="in", res=200)
ggpubr::ggarrange(fig1a, fig1b, fig1c, fig1d, ncol = 2, nrow = 2, 
                  common.legend = TRUE, labels="auto") +
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm")) 
dev.off()


### FIGURE 2
# plot biomass over distance per cruise

# scale nutrient to biomass data
coeff <- 4

fig2a <- data_figures %>%
    filter(distance > -1500) %>%
    ggplot(aes(distance, c_per_uL_mean,  col = pop, fill = pop)) + 
    geom_line(lwd = 1, position = "stack") + 
    geom_ribbon(aes(x=distance, y=c_per_uL_mean, ymin=c_per_uL_mean, ymax=c_per_uL_mean, group=pop, fill=pop), position="stack", alpha=0.5) +
    geom_point(aes(distance, NO3_NO2_mean * coeff), col = 1, pch = 16, size = 3, show.legend = FALSE) + 
    geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
    scale_fill_manual(values = pop_cols, name = "population") +
    scale_color_manual(values = pop_cols, guide = "none") +
    scale_y_continuous(name = "biomass (μgC/L)",
      sec.axis = sec_axis( trans=~./coeff, name="DIN (µmol/L)")) +
    facet_wrap( direction ~ cruise, ncol = 5) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top") +
    labs(y = "biomass (μgC/L)", x = "distance (km)")
         


png("figures/Figure_2a.png", width = 2500, height = 1500, res = 200)
print(fig2a)
dev.off()

# directional
dist_binned <- meta_gyre_d %>%
  group_by(pop, direction,
           distance = cut(distance, seq(d[1],d[2], by = res), 
                          labels = seq(d[1],d[2], by = res)[-1])) %>%
  dplyr::summarise_all(list(mean = function(x) mean(x, na.rm = TRUE), 
                            sd = function(x) sd(x, na.rm = TRUE))) %>%
  mutate(distance = as.numeric(as.character(distance)),
         pop = factor(pop, levels = c(names(pop_cols)[4], names(pop_cols)[3], names(pop_cols)[2], names(pop_cols)[1]))) %>%
  ungroup(direction)

fig2b <- dist_binned %>%
  filter(distance > -1500) %>%
  ggplot(aes(distance, c_per_uL_mean,  col = pop, fill = pop)) + 
  geom_line(lwd = 1, position = "stack") + 
  geom_ribbon(aes(x=distance, y=c_per_uL_mean, ymin=c_per_uL_mean, ymax=c_per_uL_mean, group=pop, fill=pop), position="stack", alpha=0.5) +
  geom_point(aes(distance, NO3_NO2_mean * coeff), col = 1, pch = 16, size = 3, show.legend = FALSE) + 
  geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
  scale_fill_manual(values = pop_cols, name = "population") +
  scale_color_manual(values = pop_cols, guide = "none") +
  scale_y_continuous(name = "biomass (μgC/L)",
                     sec.axis = sec_axis( trans=~./coeff, name="DIN (µmol/L)")) +
  facet_wrap(~ direction, ncol = 3) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top") +
  labs(y = "biomass (μgC/L)", x = "distance (km)")

png("figures/Figure_2b.png", width = 2500, height = 1500, res = 200)
print(fig2b)
dev.off()

png("figures/Figure_2.png", width = 2500, height = 2500, res = 200)
ggpubr::ggarrange(fig2a, fig2b, ncol = 1, nrow = 2, 
                  common.legend = TRUE, labels="auto") +
  theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm")) 
dev.off()


### FIGURE 3
# Cellular Growth

fig3 <- data_figures %>%
  filter(distance > -1500) %>%
  filter(pop == "Prochlorococcus" | pop == "Synechococcus") %>%
  drop_na(daily_growth_mean) %>%
  ggplot(aes(distance, daily_growth_mean, group = pop)) +
  geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
  geom_linerange(aes(ymin = daily_growth_mean - daily_growth_sd -  mean(daily_growth_sd, na.rm = TRUE), 
                     ymax = daily_growth_mean + daily_growth_sd +  mean(daily_growth_sd, na.rm = TRUE), col = pop), lwd= 0.5, show.legend = FALSE) +
  geom_line(aes(col= pop), lwd= 1, show.legend = TRUE) +
  scale_color_manual(values = pop_cols, name = "population") +
  facet_wrap(direction ~ cruise, ncol = 5) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top") +
  labs(y = "Net cellular growth (1/d)", x = "Distance (km)")

png("figures/Figure_3.png", width = 2500, height = 1500, res = 200)
print(fig3)
dev.off()


#directional
fig3dir <- data_figures %>%
  filter(distance > -1500) %>%
  drop_na(daily_growth_mean) %>%
  ggplot(aes(distance, daily_growth_mean, group = pop)) +
  geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
  geom_smooth(aes(col= pop), lwd= 1, show.legend = TRUE) +
  scale_color_manual(values = pop_cols, name = "population") +
  facet_wrap(~direction, ncol = 3) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top") +
  labs(y = "Net cellular growth (1/d)", x = "Distance (km)")




### FIGURE 4
# correlation plot

# directional
corr_data_north <- data_figures %>% 
  # if we want direction
  filter(direction == "north") %>%
  select(pop, c_per_uL_mean, diam_mean, daily_growth_mean, PO4_mean, NO3_NO2_mean, salinity_mean, temp_mean, daily_par_mean) %>%
  #na.omit() %>% 
  pivot_wider(names_from = pop, values_from = c(c_per_uL_mean, diam_mean, daily_growth_mean))

colnames(corr_data_north) <- c("phosphate", "nitrate", "salinity", "temperature", "daily par",
                         "biomass Pro", "biomass Syn",
                         "biomass nano", "biomass pico",
                         "diameter Pro", "diameter Syn",
                         "diameter nano", "diameter pico",
                         "growth rate Pro", "growth rate Syn",
                         "growth rate nano", "growth rate pico")

corr_data_north <- corr_data_north[,c("phosphate", "nitrate", "salinity", "temperature", "daily par",
                          "biomass Pro", "biomass Syn",
                          "biomass nano", "biomass pico",
                          "diameter Pro", "diameter Syn",
                          "growth rate Pro", "growth rate Syn")]

cor_all_north <- cor(corr_data_north, use = "complete.obs")
cor_all_p_north <- cor.mtest(corr_data_north, use = "complete.obs", conf.level = .99)

# adjust for multiple comparisons
pAdj_north <- p.adjust(c(cor_all_p_north[[1]]), method = "BH")
resAdj_north <- matrix(pAdj_north, ncol = dim(cor_all_p_north[[1]])[1])
dimnames(resAdj_north) <- dimnames(cor_all_p_north$p)

corr_data_east <- data_figures %>% 
  # if we want direction
  filter(direction == "east") %>%
  select(pop, c_per_uL_mean, diam_mean, daily_growth_mean, PO4_mean, NO3_NO2_mean, salinity_mean, temp_mean, daily_par_mean) %>%
  #na.omit() %>% 
  pivot_wider(names_from = pop, values_from = c(c_per_uL_mean, diam_mean, daily_growth_mean))

colnames(corr_data_east) <- c("phosphate", "nitrate", "salinity", "temperature", "daily par",
                               "biomass Pro", "biomass Syn",
                               "biomass nano", "biomass pico",
                               "diameter Pro", "diameter Syn",
                               "diameter nano", "diameter pico",
                               "growth rate Pro", "growth rate Syn",
                               "growth rate nano", "growth rate pico")

corr_data_east <- corr_data_east[,c("phosphate", "nitrate", "salinity", "temperature", "daily par",
                                      "biomass Pro", "biomass Syn",
                                      "biomass nano", "biomass pico",
                                      "diameter Pro", "diameter Syn",
                                      "growth rate Pro", "growth rate Syn")]

cor_all_east <- cor(corr_data_east, use = "complete.obs")
cor_all_p_east <- cor.mtest(corr_data_east, use = "complete.obs", conf.level = .99)

# adjust for multiple comparisons
pAdj_east <- p.adjust(c(cor_all_p_east[[1]]), method = "BH")
resAdj_east <- matrix(pAdj_east, ncol = dim(cor_all_p_east[[1]])[1])
dimnames(resAdj_east) <- dimnames(cor_all_p_east$p)

corr_data_south <- data_figures %>% 
  # if we want direction
  filter(direction == "south") %>%
  select(pop, c_per_uL_mean, diam_mean, daily_growth_mean, PO4_mean, NO3_NO2_mean, salinity_mean, temp_mean, daily_par_mean) %>%
  #na.omit() %>% 
  pivot_wider(names_from = pop, values_from = c(c_per_uL_mean, diam_mean, daily_growth_mean))

colnames(corr_data_south) <- c("phosphate", "nitrate", "salinity", "temperature", "daily par",
                               "biomass Pro", "biomass Syn",
                               "biomass nano", "biomass pico",
                               "diameter Pro", "diameter Syn",
                               "diameter nano", "diameter pico",
                               "growth rate Pro", "growth rate Syn",
                               "growth rate nano", "growth rate pico")

corr_data_south <- corr_data_south[,c("phosphate", "nitrate", "salinity", "temperature", "daily par",
                                      "biomass Pro", "biomass Syn",
                                      "biomass nano", "biomass pico",
                                      "diameter Pro", "diameter Syn",
                                      "growth rate Pro", "growth rate Syn")]

cor_all_south <- cor(corr_data_south, use = "complete.obs")
cor_all_p_south <- cor.mtest(corr_data_south, use = "complete.obs", conf.level = .99)

# adjust for multiple comparisons
pAdj_south <- p.adjust(c(cor_all_p_south[[1]]), method = "BH")
resAdj_south <- matrix(pAdj_south, ncol = dim(cor_all_p_south[[1]])[1])
dimnames(resAdj_south) <- dimnames(cor_all_p_south$p)

png("figures/Figure_4.png", width = 2600, height = 1200, res = 200)
par(mfrow=c(1,3))
fig_cor_north <- corrplot(cor_all_north, p.mat = resAdj_north, sig.level = 0.01, insig = "blank",
                    type = "lower", method = "color", addgrid.col = F, tl.col = "black",
                    col = colorRampPalette(c("blue", "grey90", "red"))(200))
mtext("North", at=7, line = -38, cex=1)
mtext("a", at=-1, line = -2, cex=1, font=2)
fig_cor_east <- corrplot(cor_all_east, p.mat = resAdj_east, sig.level = 0.01, insig = "blank",
                    type = "lower", method = "color", addgrid.col = F, tl.col = "black",
                    col = colorRampPalette(c("blue", "grey90", "red"))(200))
mtext("East", at=7, line = -38, cex=1)
mtext("b", at=-1, line = -2, cex=1, font=2)
fig_cor_south <- corrplot(cor_all_south, p.mat = resAdj_south, sig.level = 0.01, insig = "blank",
                    type = "lower", method = "color", addgrid.col = F, tl.col = "black",
                    col = colorRampPalette(c("blue", "grey90", "red"))(200))
mtext("South", at=7, line = -38, cex=1)
mtext("c", at=-1, line = -2, cex=1, font=2)
dev.off()









#### UNCHECKED CODE

#------------------------
# d. Supplemental Figures 
#------------------------

# corrplot of all directions
corr_data <- data_figures %>% 
  # if we want direction
  filter(direction == "south") %>%
  select(pop, c_per_uL_mean, diam_mean, daily_growth_mean, PO4_mean, NO3_NO2_mean, salinity_mean, temp_mean, daily_par_mean) %>%
  #na.omit() %>% 
  pivot_wider(names_from = pop, values_from = c(c_per_uL_mean, diam_mean, daily_growth_mean))

colnames(corr_data) <- c("phosphate", "nitrate", "salinity", "temperature", "daily par",
                         "biomass Pro", "biomass Syn",
                         "biomass nano", "biomass pico",
                         "diameter Pro", "diameter Syn",
                         "diameter nano", "diameter pico",
                         "growth rate Pro", "growth rate Syn",
                         "growth rate nano", "growth rate pico")

corr_data <- corr_data[,c("phosphate", "nitrate", "salinity", "temperature", "daily par",
                          "biomass Pro", "biomass Syn",
                          "biomass nano", "biomass pico",
                          "diameter Pro", "diameter Syn",
                          "growth rate Pro", "growth rate Syn")]

cor_all <- cor(corr_data, use = "complete.obs")
cor_all_p <- cor.mtest(corr_data, use = "complete.obs", conf.level = .99)

# adjust for multiple comparisons
pAdj <- p.adjust(c(cor_all_p[[1]]), method = "BH")
resAdj <- matrix(pAdj, ncol = dim(cor_all_p[[1]])[1])
dimnames(resAdj) <- dimnames(cor_all_p$p)

png("figures/Figure_4.png", width = 2500, height = 1600, res = 200)
fig_cor <- corrplot(cor_all, p.mat = resAdj, sig.level = 0.01, insig = "blank",
                    type = "lower", method = "color", addgrid.col = F, tl.col = "black",
                    col = colorRampPalette(c("blue", "grey90", "red"))(200))
dev.off()


# alternate form of fig2
fig2b <- data_figures %>%
  filter(distance > -1500) %>%
  ggplot(aes(distance, c_per_uL_mean,  col = pop, fill = pop)) +
  geom_line(lwd = 1, position = "fill") +
  geom_area(position = "fill",alpha = 0.5) +
  #geom_point(aes(distance, din), col = 1, pch = 1, size = 3, show.legend = FALSE) + 
  geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
  scale_fill_manual(values = pop_cols, name = "population") +
  scale_color_manual(values = pop_cols, guide = "none") +
  scale_y_continuous(name = "contribution") +
  facet_wrap(direction ~ cruise, ncol = 5) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top") +
  labs(y = "biomass (μgC/L)", x = "distance (km)")

# png("figures/Figure_2.png", width = 3500, height = 3200, res = 200)
# ggpubr::ggarrange(fig2a, fig2b, nrow = 2, common.legend = TRUE)
# dev.off()

# alternate form of fig 3
fig3 <- data_figures %>%
  filter(distance > -1500) %>%
  ggplot(aes(distance, qc_mean,  col = pop, fill = pop)) +
  geom_line(aes(group = pop), lwd = 1) +
  geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = 0.02, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
  scale_color_manual(values = pop_cols, name = "population") +
  scale_y_continuous(trans='log10') +
  facet_wrap(. ~ cruise) +
  theme_bw(base_size = 20) +
  labs(y = "carbon quota (pg C/cell)", x = "distance (km)")

# nutrients
data_nutr <- data_figures[, c("cruise", "pop", "distance", "NO3_NO2_mean", "PO4_mean")]
data_nutr <- data_nutr %>%
  dplyr::rename(NO3_NO2 = NO3_NO2_mean, PO4 = PO4_mean) %>%
  gather(nutrient, concentration, 4:5)
data_nutr$modeled <- "modeled"
actual <- which(data_nutr$cruise == "KOK1606" | 
                  data_nutr$cruise == "MGL1704" | 
                  data_nutr$cruise == "KM1906" | 
                  data_nutr$cruise == "TN398")
data_nutr$modeled[actual] <- "observed"
ind_nutr <- which(data_nutr$cruise == "KM1923" & data_nutr$distance > 1750)
data_nutr <- data_nutr[-ind_nutr,]

nutrient_names <- list("NO3_NO2" = "nitrate", "PO4" = "phosphate")
nutrient_labeller <- function(variable, value){
  return(nutrient_names[value])
}


fig_nutr <- data_nutr %>%
  ggplot(aes(distance, concentration, col = cruise, shape = modeled)) + 
  geom_point(size = 2) +
  scale_color_manual(values=viridis(9, alpha = 1, begin = 0, end = 1, direction = 1)) +
  theme_bw(base_size = 15) +
  facet_wrap(. ~ nutrient, scale = "free", labeller=nutrient_labeller) +
  labs(y = "nutrient concentration (μg / L)", x = "distance (km)")

png(paste0("figures/", "nutr-cruise.png"), width = 2500, height = 1200, res = 200)
print(fig_nutr)
dev.off()



### plotting nutrients against phytoplankton diameter
din <- data_figures %>%
  select(cruise, pop, distance, contains('NO3')) %>%
  rename(mean = contains("mean")) %>%
  filter(!is.na(mean)) %>%
  ggplot(aes(distance, mean,  col = cruise)) + 
  geom_rect(aes(xmin = mean(front_uncertainties$down), xmax = mean(front_uncertainties$up), ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey") +
  geom_point(aes(group = cruise)) +
  #geom_line(aes(group = cruise), lwd = 0.5) +
  theme_bw(base_size = 20) +
  labs(y = "Nitrate and Nitrite (μmol/L)", x = "distance (km)")


### Temp, nitrate, PAR
temp <- data_figures %>%
  select(cruise, pop, distance, contains('temp')) %>%
  rename(mean = contains("mean")) %>%
  filter(!is.na(mean)) %>%
  ggplot(aes(distance, mean,  col = cruise)) + 
  geom_rect(aes(xmin = mean(front_uncertainties$down), xmax = mean(front_uncertainties$up), ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey") +
  geom_point(aes(group = cruise)) +
  geom_line(aes(group = cruise), lwd = 0.5) +
  theme_bw(base_size = 20) +
  labs(y = "Temperature (ºC)", x = "distance (km)")

par <- data_figures %>%
  select(cruise, pop, distance, contains('daily_par')) %>%
  rename(mean = contains("mean")) %>%
  filter(!is.na(mean)) %>%
  ggplot(aes(distance, mean,  col = cruise)) + 
  geom_rect(aes(xmin = mean(front_uncertainties$down), xmax = mean(front_uncertainties$up), ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey") +
  geom_point(aes(group = cruise)) +
  geom_line(aes(group = cruise), lwd = 0.5) +
  theme_bw(base_size = 20) +
  labs(y = "daily PAR (µE m-2)", x = "distance (km)")


png(paste0("figures/Figure_3.png"), width = 3500, height = 1600, res = 200)
ggpubr::ggarrange(temp, par, din, nrow = 1, common.legend = TRUE)
dev.off()


### growth rate

fig4 <- data_figures %>%
  filter(distance > -1500) %>%
  filter(!is.na(daily_growth_mean)) %>%
  ggplot(aes(distance, daily_growth_mean,  col = pop, fill = pop)) +
  geom_line(aes(group = pop), lwd = 1) +
  geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
  scale_color_manual(values = pop_cols, name = "population") +
  facet_wrap(. ~ cruise) +
  theme_bw(base_size = 20) +
  labs(y = "growth rate", x = "distance (km)")

png(paste0("figures/","Figure_4.png"), width = 2500, height = 2000, res = 200)
print(fig4)
dev.off()

### MLD
data_figures$Month <- lubridate::month(data_figures$date_mean)

fig_mld <- data_figures %>%
  ggplot(aes(distance, MLD_mean,  col = as.factor(Month))) +
  geom_point(aes(group = pop), size = 2) +
  geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
  facet_wrap(. ~ cruise) +
  theme_bw(base_size = 20) +
  labs(y = "mixed layer depth (m)", x = "distance (km)") +
  scale_y_reverse() +
  scale_color_discrete()

png(paste0("figures/","mld_seasonally.png"), width = 2500, height = 2000, res = 200)
print(fig_mld)
dev.off()

### MLD and NO3
data_figures$MLD_NO3 <- data_figures$MLD_mean * data_figures$NO3_NO2_mean

fig_mld_no3 <- data_figures %>%
  ggplot(aes(distance, MLD_NO3)) +
  geom_point(aes(group = pop), size = 2) +
  geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
  facet_wrap(. ~ cruise) +
  theme_bw(base_size = 20) +
  labs(y = "MLD * DIN", x = "distance (km)")

png(paste0("figures/","mld_no3.png"), width = 2500, height = 2000, res = 200)
print(fig_mld_no3)
dev.off()

