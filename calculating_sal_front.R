################################
Salinity Front Code
Jordan Winter
Last edited 5/6/22
################################


########## loading packages and data
# load packages

library(tidyverse)
library(maps)
library(geosphere)

# load data
data <- read.csv("~/Downloads/paper/all_data_PROGRESS.csv")
data$hour <- as.POSIXct(data$hour, format = "%Y-%m-%d %H", tz = "GMT")

########### getting cruise names
cruises <- unique(data$cruise)

########### creating a loop to get the second derivative of salinity and finding where it is 0

# create empty dataframe
sal_front <- NULL

# loop through each cruise
for (cruise_name in cruises){
	
	######### subset the data for one cruise
	data_small <- subset(data, cruise == cruise_name)
	
	# make sure data is in the right time order and get hourly means
	data_small <- data_small[order(data_small$hour),]
	hourly <- group_by(data_small, hour)
	data_small <- dplyr::summarize_all(hourly, mean)
	data_small <- as.data.frame(data_small)
	data_small$hour <- as.POSIXct(data_small$hour, format = "%Y-%m-%d %H", tz = "GMT")
	
	# make sure there is a row for every hour in the dataframe
	data_small <- data_small %>% complete(hour = seq.POSIXt(min(hour), max(hour), by = "hour")) # get all hours
	data_small <- as.data.frame(data_small)
	
	# get rid of NaNs
	data_small <- subset(data_small, is.na(lat) == FALSE)
	
	######## calculate the first derivative
	deriv1 <- NULL
	deriv1$first <- diff(data_small$salinity)
	deriv1 <- as.data.frame(deriv1)
	i <- nrow(data_small) - 1
	deriv1$time <- tail(data_small$hour, n = i)
	deriv1$time <- as.POSIXct(deriv1$time, format = "%Y-%m-%d %H", tz = "GMT")
	
	# get rid of outliers with large improbable changes in salinity
	deriv1 <- subset(deriv1, first < 2 & first > -2)
	
	# make sure there is a row for every hour in the dataframe
	deriv1 <- deriv1 %>% complete(time = seq.POSIXt(min(time), max(time), by = "hour"))
	
	######### calculate the second derivative
	deriv2 <- NULL
	deriv2$second <- diff(deriv1$first)
	deriv2 <- as.data.frame(deriv2)
	i <- nrow(deriv1) - 1
	deriv2$time <- tail(deriv1$time, n = i)
	deriv2$time <- as.POSIXct(deriv2$time, format = "%Y-%m-%d %H", tz = "GMT")
	
	# take the absolute value so all values are positive and it is easier to find the minimum
	deriv2$second <- abs(deriv2$second)
	
	########## find the value closest to 0 for the second derivative
	min <- min(deriv2$second, na.rm = T)
	ind <- which(deriv2$second == min)
	
	# take the time of the first index with the value closest to 0 as the salinity front
	time_front <- deriv2$time[ind[1]]
	time_front <- as.character(time_front)
	
	########### for cruises that found a different salinity front, use these times as the front
	if (cruise_name == "Gradients 4"){
		time_front <- "2021-12-10 11:00:00"
	}
	if (cruise_name == "TN398"){
		time_front <- "2021-12-27 02:00:00"
	}
	if (cruise_name == "SR1917"){
		time_front <- "2019-11-25 06:00:00"
	}
	if (cruise_name == "KM1923"){
		time_front <- "2019-12-12 09:00:00"
	}
	if (cruise_name == "KM1713"){
		time_front <- "2017-09-22 03:00:00"
	}
	if (cruise_name == "Gradients 1"){
		time_front <- "2016-04-25 09:00:00"
	}
	if (cruise_name == "Gradients 3"){
		time_front <- "2019-04-21 23:00:00"
	}
	
	########### get the lat and lon of the salinity front
	data_small$hour <- as.character(data_small$hour)
	ind <- which(data_small$hour == time_front)
	lat_zero <- data_small$lat[ind]
	lon_zero <- data_small$lon[ind]
	
	########## define inside vs outside the gyre
	data_small$gyre <- "not"
	
	# take a random point in the gyre as a metric for where the salinity front is in space in relation to the gyre
	gyre_lat <- 22.75
	gyre_lon <- -160
	
	# if the lat and lon are greater than the gyre lat and lon, mark specific regions as gyre
	if (lat_zero > gyre_lat & lon_zero > gyre_lon){
		ind <- which(data_small$lat < lat_zero & lon_zero > gyre_lon)
		data_small$gyre[ind] <- "gyre"
		ind <- which(data_small$lat == lat_zero)
		data_small$gyre[ind] <- "zero"
	}
	
	# if the lat is less than the gyre lat and the lon is greater than the gyre lon, mark certain regions as gyre
	if (lat_zero < gyre_lat & lon_zero > gyre_lon){
		ind <- which(data_small$lat > lat_zero & lon_zero > gyre_lon)
		data_small$gyre[ind] <- "gyre"
		ind <- which(data_small$lat == lat_zero)
		data_small$gyre[ind] <- "zero"
	}
	
	# add in the second and third salinity fronts of Gradients 4
	if (cruise_name == "Gradients 4"){
		time_front <- "2021-11-20 07:00:00"
		ind <- which(data_small$hour == time_front)
		lat_zero <- data_small$lat[ind]
		lon_zero <- data_small$lon[ind]
		data_small$gyre[ind] <- "zero"
		ind <- which(data_small$lat > lat_zero & lon_zero > gyre_lon)
		data_small$gyre[ind] <- "not"
		time_front <- "2021-11-30 13:00:00"
		ind <- which(data_small$hour == time_front)
		lat_zero <- data_small$lat[ind]
		lon_zero <- data_small$lon[ind]
		data_small$gyre[ind] <- "zero"
	}
	
	# making the SR1917 cruise only have one front
	if (cruise_name == "SR1917"){
		ind <- which(data_small$lon > 0 | data_small$lon < -150)
		data_small$gyre[ind] <- "gyre"
	}
	data_small$cruise <- cruise_name
	
	########## add to the larger dataframe
	sal_front <- rbind(sal_front, data_small)
}

########## create map of gyre boundary to see if it looks reasonable
worldmap <- map_data("world")
worldmap <- subset(worldmap, long > -225 & long < -110)
worldmap <- subset(worldmap, lat > -10 & lat < 70)

ind <- which(sal_front$lon > 0)
sal_front$lon[ind] <- sal_front$lon[ind] - 360
sal_front <- subset(sal_front, lon > -250)

png(paste0("~/Downloads/paper/figures/sal_map.png"),width=18, height=10, unit="in", res=200)
g <- ggplot() +
	 geom_map(data = worldmap, map = worldmap, aes(long, lat, map_id = region), color = "white", fill = "lightgray") +
	 geom_point(data = sal_front, aes(lon, lat, color = gyre), size=2, alpha = 0.7, show.legend = T) +
	 labs(x='Longitude', y= 'Latitude') +
	 theme_bw() +
	 scale_color_manual(values = viridis(3)) +
	 theme(text = element_text(size = 15)) +
	 ggtitle("Cruise Tracks vs Salinity Front")
print(g)
dev.off()

############## making the gyre boundary a 0 point and finding the distance from there

# add a column for the distance to the dataframe
sal_front$haversine <- NA

# for cruises with one salinity front, loop through each cruise and calculate the distance using distHaversine()
sal_fixed <- NULL

# make sure cruise list only has cruises with one salinity front and negative lon values
cruises <- cruises[cruises != "Gradients 4" & cruises != "SR1917"]

# loop through each cruise
for (cruise_name in cruises){
	
	print(cruise_name)
	
	# subset data only fron one cruise
	sal_small <- subset(sal_front, cruise == cruise_name)
	
	# initialize dataframe
	D <- NULL
	
	# make sure lat and lon are values
	sal_small <- subset(sal_small, is.na(lat) == F)
	sal_small <- subset(sal_small, is.na(lon) == F)
	
	# find the index in the dataframe for the salinity front
	zero <- which(sal_small$gyre == "zero")
	
	# calculate the distance from the salinity front to the lat and lon of each row of the dataframe
	for (i in 1:nrow(sal_small)){
		d <- distHaversine(sal_small[c(zero[1],i), 4:3])
		sal_small$haversine[i] <- d 		#distance is in meters
	}
	
	# if the location is not in the gyre, make sure the distance is positive
	for (i in 1:nrow(sal_small)){
		if (sal_small$gyre[i] == "not"){
			sal_small$haversine[i] <- abs(sal_small$haversine[i])
		}
	
	# if the location is in the gyre, make sure the distance is negative
		if (sal_small$gyre[i] == "gyre"){
			sal_small$haversine[i] <- abs(sal_small$haversine[i]) * -1
		}
	}
	
	# add all this information to the larger dataframe
	sal_fixed <- rbind(sal_fixed, sal_small)
}

# only save the rows needed
sal_fixed <- sal_fixed[c("hour", "lat", "lon", "cruise", "gyre", "haversine")]

########### finding the distance using time instead of lat and lon

# for Gradients 4 (multiple salinity fronts) and SR1917 (positive lon values)
cruises <- c("Gradients 4", "SR1917")

# loop through each cruise
for (cruise_name in cruises){
	
	print(cruise_name)
	
	# subset larger dataframe
	sal_small <- subset(sal_front, cruise == cruise_name)
	
	# make the lon values positive if they're less than -180
	test <- which(sal_small$lon < -180)
	if (length(test) != 0){
		sal_small$lon[test] <- sal_small$lon[test] + 360
	}
	
	D <- NULL
	
	# make sure the dataframe only consists of rows that have lat and lon
	sal_small <- subset(sal_small, is.na(lat) == F)
	sal_small <- subset(sal_small, is.na(lon) == F)
	
	# calculate the distance between each neighboring lat and lon and save it
	for (i in 2:nrow(sal_small)){
		d <- distHaversine(sal_small[(i-1):i, 4:3])
		D <- c(D, d)
		flush.console()
	}
	
	# make sure the distance dataframe is the same length as the subsetted dataframe
	D[length(D) + 1] <- 0
	
	# add the distances up along the cruise track
	sal_small$dist <- D
	sal_small$dist_add <- NULL
	sal_small$dist_add[1] <- sal_small$dist[1]
	for (i in 2:nrow(sal_small)){
		sal_small$dist_add[i] <- sal_small$dist_add[i-1] + sal_small$dist[i]
	}
	
	# remove the added 0 value that make the dataframe the same length as the subsetted dataframe
	sal_small <- subset(sal_small, dist != 0)
	
	# find where the salinity front is
	zero_km <- which(sal_small$gyre == "zero")
	
	# if there is one salinity front, subtract the distance at the salinity front from all other distances
	if (length(zero_km == 1)){
	zero <- sal_small$dist_add[zero_km]
	sal_small$haversine <- zero - sal_small$dist_add
	}
	
	# if there are multiple salinity fronts, find the nearest salinity front to subtract the distance from
	if (length(zero_km > 1)){
		for (i in 1:nrow(sal_small)){
			zero <- sal_small$dist_add[zero_km]
			sal_small$haversine[i] <- min(abs(zero - sal_small$dist_add[i]))
		}
	}
	
	# make sure values not in the gyre are positive and values in the gyre are negative
	for (i in 1:nrow(sal_small)){
		if (sal_small$gyre[i] == "not"){
			sal_small$haversine[i] <- abs(sal_small$haversine[i])
		}
		if (sal_small$gyre[i] == "gyre"){
			sal_small$haversine[i] <- abs(sal_small$haversine[i]) * -1
		}
	}
	
	# only add the columns needed to the dataframe
	sal_small <- sal_small[c("hour", "lat", "lon", "cruise", "gyre", "haversine")]
	sal_fixed <- rbind(sal_fixed, sal_small)
}

# make a csv of the data
write.csv(sal_fixed, "~/Downloads/paper/sal_front.csv")
