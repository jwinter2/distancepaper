library(tidyverse)
library(arrow)
library(geosphere)
library(sf)
library(plotly)
library(viridis)
library(dplyr)

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

# change Gradients cruise names to their actual cruise names
id_grad1 <- which(meta$cruise == "Gradients 1")
meta$cruise[id_grad1] <- "KOK1606"
id_grad2 <- which(meta$cruise == "Gradients 2")
meta$cruise[id_grad2] <- "MGL1704"
id_grad3 <- which(meta$cruise == "Gradients 3")
meta$cruise[id_grad3] <- "KM1906"
id_grad4 <- which(meta$cruise == "Gradients 4")
meta$cruise[id_grad4] <- "TN397"

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
    dplyr::group_by(cruise) %>% 
    dplyr::summarize(sal_front = max(sal_front, na.rm = TRUE),
      sal_front_up = max(sal_front_up, na.rm = TRUE),
      sal_front_down = max(sal_front_down, na.rm = TRUE)) 

 
### plot locations of all abrupt changes of salinity
meta_distance_binned %>%  ggplot() +
  geom_point(aes(distance, salinity)) +
  geom_point(aes(distance, smooth_salinity), col = "red") +
  geom_point(aes(distance, smooth_salinity_down), col = "green") +
  geom_point(aes(distance, smooth_salinity_up), col = "blue") +
  geom_hline(data = tz, aes( yintercept = sal_front), col = "red") +
  geom_hline(data = tz, aes( yintercept = sal_front_down), col = "green") +
  geom_hline(data = tz, aes( yintercept = sal_front_up), col = "blue") +
  facet_wrap(. ~ cruise, scale = "free_x")


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
# the date is to make sure this is not finding a difference between different cruises
id <- which(diff(meta_gyre$gyre) != 0 & diff(meta_gyre$date) < 345600)
boundaries <- meta_gyre[id,] %>% 
  dplyr::group_by(date) %>% 
  dplyr::summarize(lat = mean(lat), 
    lon = mean(lon),
    cruise = unique(cruise))

id_down <- which(diff(meta_gyre$gyre_down) != 0 & diff(meta_gyre$date) < 345600)
boundaries_down <- meta_gyre[id_down,] %>% 
  dplyr::group_by(date) %>% 
  dplyr::summarize(lat = mean(lat), 
    lon = mean(lon),
    cruise = unique(cruise))

id_up <- which(diff(meta_gyre$gyre_up) != 0 & diff(meta_gyre$date) < 345600)
boundaries_up <- meta_gyre[id_up,] %>% 
  dplyr::group_by(date) %>% 
  dplyr::summarize(lat = mean(lat), 
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

# plot cruise tracks
g <- plot_geo(meta_gyre_d, lat = ~ lat, lon = ~ lon, color = ~ cruise, mode = "scatter", colors = viridis(9, alpha = 1, begin = 0, end = 1, direction = 1)) %>%
  layout(geo = geo)
#g


#--------------------------
# b. PLOTTING OVER DISTANCE
#--------------------------
### Calculate mean and sd over binned distance from the edges of the NPSG
res <- 100 # km (i.e data binned every 100 km)


### remove data far outside the gyre
faraway_id <- which(meta_gyre_d$distance > 2000)
meta_gyre_d <- meta_gyre_d[-faraway_id,]

d <- range(meta_gyre_d$distance)
data_figure <- meta_gyre_d %>%
  group_by(cruise, pop, 
    distance = cut(distance, seq(d[1],d[2], by = res), 
      labels = seq(d[1],d[2], by = res)[-1])) %>%
  dplyr::summarise_all(list(mean = function(x) mean(x, na.rm = TRUE), 
    sd = function(x) sd(x, na.rm = TRUE))) %>%
mutate(distance = as.numeric(as.character(distance))) %>%
ungroup(cruise)

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

for (cruise_name in unique(front_uncertainties$cruise)){
  small_unc <- front_uncertainties %>% filter(cruise == cruise_name)
  up <- which(meta_gyre_d$cruise == cruise_name & meta_gyre_d$distance <= small_unc$up
              & meta_gyre_d$distance >= 0)
  meta_gyre_d$gyre[up] <- "transition"
  down <- which(meta_gyre_d$cruise == cruise_name & meta_gyre_d$distance >= small_unc$down
              & meta_gyre_d$distance <= 0)
  meta_gyre_d$gyre[down] <- "transition"
}

#plots on flat map projection
for_map <- meta_gyre_d
for_map$lon <- for_map$lon - 360
worldmap <- map_data("world")
worldmap <- subset(worldmap, long > -170 & long < -110)
worldmap <- subset(worldmap, lat > -10 & lat < 70)

#choose a parameter
type <- "cruise"
type <- "gyre"
type <- "salinity"

png(paste0("figures/", type, "_map.png"),width=12, height=10, unit="in", res=200)
g <- ggplot() +
  geom_map(data = worldmap, map = worldmap, aes(long, lat, map_id = region), color = "white", fill = "lightgray") +
  geom_point(data = for_map, aes(lon, lat, color = gyre), size=2, alpha = 0.7, show.legend = T) +
  theme_bw() +
  #scale_color_viridis() +
  scale_color_manual(values = viridis(3)) +
  theme(text = element_text(size = 20))
print(g)
dev.off()


### plotting nutrients
data_nutr <- data_figure[, c("cruise", "pop", "distance", "NO3_NO2_mean", "PO4_mean")]
data_nutr <- data_nutr %>%
  dplyr::rename(NO3_NO2 = NO3_NO2_mean, PO4 = PO4_mean) %>%
  gather(nutrient, concentration, 4:5)
data_nutr$modeled <- "modeled"
actual <- which(data_nutr$cruise == "KOK1606" | data_nutr$cruise == "MGL1704" | data_nutr$cruise == "KM1906" | data_nutr$cruise == "TN398")
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

### choose a parameter
para <- "n_per_uL"; ylab <- "abundance (cells / μL)"; name <- "abundance"
# or
para <- "c_per_uL"; ylab <- "biomass (μgC / L)"; name <- "biomass"
# or
para <- "diam"; ylab <- "equivalent spherical diameter (μm)"; name <- "diameter"


### set colors
pop_cols <- c("prochloro" = rocket(7)[6], "synecho" = rocket(7)[4], "picoeuk" = rocket(7)[2])

### pop as factors
data_figure$pop <- factor(data_figure$pop, levels = c("prochloro", "synecho", "picoeuk", "croco", "beads", "unknown"))

### plotting parameter over distance per cruise
fig3 <- data_figure %>%
 select(cruise, pop, distance, contains(para)) %>%
 rename(mean = contains("mean")) %>%
  ggplot(aes(distance, mean,  col = pop, fill = pop)) + 
    geom_line(aes(group = pop), lwd = 1, position = "stack") + 
    geom_area(position = "stack", stat = "identity", alpha = 0.5) +
    geom_rect(data = front_uncertainties, aes(xmin = down, xmax = up, ymin = -Inf, ymax = Inf), alpha= 0.25, inherit.aes = FALSE) +
    scale_fill_manual(values = pop_cols, name = "population", labels = c("prochlorococus", "synechococcus", "picoeukaryotes")) +
    scale_color_manual(values = pop_cols, guide = "none") +
    facet_wrap(. ~ cruise, scale = "free_x") +
    theme_bw(base_size = 20) +
    labs(y = ylab, x = "distance (km)")


### save plot
png(paste0("figures/",name,"-distance-stacked.png"), width = 2500, height = 1600, res = 200)
print(fig3)
dev.off()


### plotting nutrients against phytoplankton diameter

### take mean diameter per day
meta_gyre_d$day <- substr(meta_gyre_d$date, 1, 10)
diam_data <- meta_gyre_d %>%
  group_by(cruise, pop, day, gyre) %>%
  dplyr::summarise_at(c("diam", "c_per_uL", "NO3_NO2"), list(mean), na.rm=T)

diam_data$pop <- factor(diam_data$pop, levels = c("prochloro", "synecho", "picoeuk", "croco", "beads", "unknown"))
pop_names <- list("prochloro" = "prochlorococcus", "synecho" = "synechococcus", "picoeuk" = "picoeukaryotes")
pop_labeller <- function(variable, value){
    return(pop_names[value])
}

fig4 <- diam_data %>%
  ggplot(aes(diam, NO3_NO2, col = gyre, levels = pop)) + 
  geom_point() +
  theme_bw(base_size = 20) +
  scale_color_manual(values=c(viridis(5)[1], viridis(5)[4], viridis(5)[5])) +
  facet_grid(.~ pop, scale = "free", labeller = pop_labeller) +
  labs(x = "equivalent spherical diameter (μm)", y = "nutrient concentration (μg/L)")


### save plot

png(paste0("figures/","nutr-diam.png"), width = 2500, height = 800, res = 200)
print(fig4)
dev.off()

### correlation plot
corplot_all_df <- meta_gyre_d[, c("pop", "c_per_uL", "diam", "NO3_NO2", "salinity", "temp", "par")]
corplot_all_df <- corplot_all_df %>% na.omit() %>% pivot_wider(names_from = pop, values_from = c(c_per_uL, diam))
colnames(corplot_all_df) <- c("nitrate", "salinity", "temperature", "par", "biomass pico", "biomass syn", "biomass pro", "diameter pico", "diameter syn", "diameter pro")

cor_all <- cor(corplot_all_df, use = "complete.obs")
cor_all_p <- cor.mtest(corplot_all_df, use = "complete.obs", conf.level = .99)

png(paste0("figures/","all-corr.png"), width = 2500, height = 1600, res = 200)
fig_cor <- corrplot(cor_all, p.mat = cor_all_p$p, sig.level = 0.01, insig = "blank", type = "lower", method = "square", addgrid.col = F, tl.col = "black", col = c('steelblue3', 'firebrick3'))
dev.off()
