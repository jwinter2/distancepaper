################################
Data Retrieval and Cleaning Code
Jordan Winter
Last edited 4/26/22
################################

######### load packages

library(tidyverse)
library(plyr)
library(arrow)


########## function for reading sfl file
read_sfl <- function(x){
  df <- read_delim(x, delim="\t", guess_max = 2000) # for case when the first 2000 rows are NA's

    #parse cruise name and serial number of instrument
    exp <- unlist(list(strsplit(sub(".sfl", "", basename(x)),"_")))

      if(length(exp) > 2) { cruise <- paste(exp[1],exp[2],sep="_")
      } else if(length(exp) ==2) cruise <- exp[1]
      print(cruise)
      inst <-  sub(".sfl","",exp[length(exp)])

  return(df)
}


########## read data

# TN398 cruise SeaFlow data
tn398_sfl <- read.csv("~/Downloads/paper/raw_data/TN398_SeaFlow_stat_table.csv")
tn398_sfl$time <- as.POSIXct(tn398_sfl$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")

# TN398 temp, salinity, PAR data
tn398_env <- read_sfl("~/Downloads/paper/raw_data/TN398_740.sfl")
tn398_env <- as.data.frame(tn398_env)
tn398_env$time <- as.POSIXct(tn398_env$DATE, format = "%Y-%m-%d %H:%M", tz = "GMT")
tn398_env <- tn398_env[, c("LAT", "LON", "CONDUCTIVITY", "SALINITY", "OCEAN TEMP", "PAR", "time")]

# TN398 nutrient data from the Marine Chemistry Lab at UW
tn398_nutr <- read.csv("~/Downloads/paper/raw_data/TN398_nutr.csv")

# Gradients 1, 2, 3, 4 SeaFlow data from CMAP
grad1_sfl <- read.csv("~/Downloads/paper/raw_data/SCOPE_16_SeaFlow_stat_table.csv")
grad1_sfl$time <- as.POSIXct(grad1_sfl$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
grad1_sfl$cruise <- "Gradients 1"

grad2_sfl <- read.csv("~/Downloads/paper/raw_data/MGL1704_SeaFlow_stat_table.csv")
grad2_sfl$time <- as.POSIXct(grad2_sfl$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
grad2_sfl$cruise <- "Gradients 2"

grad3_sfl <- read.csv("~/Downloads/paper/raw_data/KM1906_SeaFlow_stat_table.csv")
grad3_sfl$time <- as.POSIXct(grad3_sfl$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
grad3_sfl$cruise <- "Gradients 3"

grad4_sfl <- read.csv("~/Downloads/paper/raw_data/TN397_740_SeaFlow_stat_table.csv")
grad4_sfl$time <- as.POSIXct(grad4_sfl$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
grad4_sfl$cruise <- "Gradients 4"

grad_sfl <- rbind(grad1_sfl, grad2_sfl)
grad_sfl <- rbind(grad_sfl, grad3_sfl)
grad_sfl <- rbind(grad_sfl, grad4_sfl)

# Gradients 1, 2, 3 temp, salinity, PAR data from parquet files
grad1_env <- arrow::read_parquet("~/Downloads/sfl_stats/KOK1606/SeaFlow_KOK1606_v1.3_2020-08-21.parquet")
grad1_env <- as.data.frame(grad1_env)
grad1_env <-  subset(grad1_env, flag == 0)
grad1_env <- subset(grad1_env, quantile == 2.5)
grad1_env$time <- as.POSIXct(grad1_env$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
grad1_env$cruise <- "Gradients 1"

grad2_env <- read_parquet("~/Downloads/sfl_stats/MGL1704/SeaFlow_MGL1704_v1.3_2020-08-21.parquet")
grad2_env <- as.data.frame(grad2_env)
grad2_env <-  subset(grad2_env, flag == 0)
grad2_env <- subset(grad2_env, quantile == 2.5)
grad2_env$time <- as.POSIXct(grad2_env$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
grad2_env$cruise <- "Gradients 2"

grad3_env <- read_parquet("~/Downloads/sfl_stats/KM1906/SeaFlow_KM1906_v1.3_2020-08-21.parquet")
grad3_env <- as.data.frame(grad3_env)
grad3_env <-  subset(grad3_env, flag == 0)
grad3_env <- subset(grad3_env, quantile == 2.5)
grad3_env$time <- as.POSIXct(grad3_env$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
grad3_env$cruise <- "Gradients 3"

grad4_env <- read_sfl("~/Downloads/curation/TN397_740/TN397_740.sfl")
grad4_env <- as.data.frame(grad4_env)
grad4_env <-  subset(grad4_env, is.na(LAT) == F)
grad4_env$cruise <- "Gradients 4"
grad4_env <- grad4_env[,c("DATE", "LAT", "LON", "PAR", "SALINITY", "OCEAN TEMP", "cruise")]
colnames(grad4_env) <- c("time", "lat", "lon", "par", "salinity", "temp", "cruise")

grad_env <- rbind(grad1_env, grad2_env)
grad_env <- rbind(grad_env, grad3_env)
grad_env <- grad_env[,c("time", "lat", "lon", "par", "salinity", "temp", "cruise")]
grad_env <- rbind(grad_env, grad4_env)

# Gradients 1, 2, 3 nutrient data from CMAP
# https://simonscmap.com/catalog/datasets/KOK1606_Gradients1_Nutrients
# https://simonscmap.com/catalog/datasets/MGL1704_Gradients2_Nutrients
grad_1_nutr <- read.csv("~/Downloads/paper/raw_data/grad1_nutr.csv")
grad_1_nutr$time <- as.POSIXct(grad_1_nutr$time, format = "%Y-%m-%d %H:%M", tz = "GMT")

grad_2_nutr <- read.csv("~/Downloads/paper/raw_data/grad2_nutr.csv")
grad_2_nutr$time <- as.POSIXct(grad_2_nutr$time, format = "%Y-%m-%d %H:%M", tz = "GMT")

grad_3_nutr <- read.csv("~/Downloads/paper/raw_data/grad3_nutr.csv")
grad_3_nutr$time <- as.POSIXct(grad_3_nutr$time, format = "%Y-%m-%d %H:%M", tz = "GMT")

grad_nutr <- rbind(grad_1_nutr, grad_2_nutr)
grad_nutr <- rbind(grad_nutr, grad_3_nutr)

# KM1712, KM1713, KM1923, KM1502, SR1917, TN271 SeaFlow data from CMAP
km1712_sfl <- read.csv("~/Downloads/paper/raw_data/KM1712_SeaFlow_stat_table.csv")
km1712_sfl$time <- as.POSIXct(km1712_sfl$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
km1712_sfl$cruise <- "KM1712"

km1713_sfl <- read.csv("~/Downloads/paper/raw_data/KM1713_SeaFlow_stat_table.csv")
km1713_sfl$time <- as.POSIXct(km1713_sfl$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
km1713_sfl$cruise <- "KM1713"

km1502_sfl <- read.csv("~/Downloads/paper/raw_data/SCOPE_2_SeaFlow_stat_table.csv")
km1502_sfl$time <- as.POSIXct(km1502_sfl$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
km1502_sfl$cruise <- "KM1502"

km1923_sfl <- read.csv("~/Downloads/paper/raw_data/KM1923_751_SeaFlow_stat_table.csv")
km1923_sfl$time <- as.POSIXct(km1923_sfl$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
km1923_sfl$cruise <- "KM1923"

sr1917_sfl <- read.csv("~/Downloads/paper/raw_data/SR1917_SeaFlow_stat_table.csv")
sr1917_sfl$time <- as.POSIXct(sr1917_sfl$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
sr1917_sfl$cruise <- "SR1917"

tn271_sfl <- read.csv("~/Downloads/paper/raw_data/Thompson_9_SeaFlow_stat_table.csv")
tn271_sfl$time <- as.POSIXct(tn271_sfl$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
tn271_sfl$cruise <- "TN271"

other_sfl <- rbind(km1712_sfl, km1713_sfl)
other_sfl <- rbind(other_sfl, km1502_sfl)
other_sfl <- rbind(other_sfl, km1923_sfl)
other_sfl <- rbind(other_sfl, sr1917_sfl)
other_sfl <- rbind(other_sfl, tn271_sfl)
other_sfl <- rbind(other_sfl, km1923_sfl)

# KM1712, KM1713, KM1923, KM1502, SR1917, TN271 temp, salinity, PAR data from parquet files
km1712_env <- arrow::read_parquet("~/Downloads/sfl_stats/KM1712/SeaFlow_KM1712_v1.3_2020-08-21.parquet")
km1712_env <- as.data.frame(km1712_env)
km1712_env <-  subset(km1712_env, flag == 0)
km1712_env <- subset(km1712_env, quantile == 2.5)
km1712_env$time <- as.POSIXct(km1712_env$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
km1712_env$cruise <- "KM1712"

km1713_env <- arrow::read_parquet("~/Downloads/sfl_stats/KM1713/SeaFlow_KM1713_v1.3_2020-08-21.parquet")
km1713_env <- as.data.frame(km1713_env)
km1713_env <-  subset(km1713_env, flag == 0)
km1713_env <- subset(km1713_env, quantile == 2.5)
km1713_env$time <- as.POSIXct(km1713_env$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
km1713_env$cruise <- "KM1713"

sr1917_env <- arrow::read_parquet("~/Downloads/sfl_stats/SR1917/SeaFlow_SR1917_v1.3_2020-08-21.parquet")
sr1917_env <- as.data.frame(sr1917_env)
sr1917_env <-  subset(sr1917_env, flag == 0)
sr1917_env <- subset(sr1917_env, quantile == 2.5)
sr1917_env$time <- as.POSIXct(sr1917_env$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
sr1917_env$cruise <- "SR1917"

km1502_env <- arrow::read_parquet("~/Downloads/sfl_stats/KM1502/SeaFlow_KM1502_v1.3_2020-08-21.parquet")
km1502_env <- as.data.frame(km1502_env)
km1502_env <-  subset(km1502_env, flag == 0)
km1502_env <- subset(km1502_env, quantile == 2.5)
km1502_env$time <- as.POSIXct(km1502_env$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
km1502_env$cruise <- "KM1502"

tn271_env <- arrow::read_parquet("~/Downloads/sfl_stats/TN271/SeaFlow_TN271_v1.3_2020-08-21.parquet")
tn271_env <- as.data.frame(tn271_env)
tn271_env <-  subset(tn271_env, flag == 0)
tn271_env <- subset(tn271_env, quantile == 2.5)
tn271_env$time <- as.POSIXct(tn271_env$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
tn271_env$cruise <- "TN271"

km1923_env <- read_sfl("~/Downloads/paper/raw_data/KM1923_751_751.sfl")
km1923_env <- as.data.frame(km1923_env)
km1923_env$time <- as.POSIXct(km1923_env$DATE, format = "%Y-%m-%d %H:%M", tz = "GMT")
km1923_env <- km1923_env[, c("LAT", "LON", "SALINITY", "OCEAN TEMP", "PAR", "DATE")]
colnames(km1923_env) <- c("lat", "lon", "salinity", "temp", "par", "time")
km1923_env$cruise <- "KM1923"

other_env <- rbind(km1712_env, km1713_env)
other_env <- rbind(other_env, sr1917_env)
other_env <- rbind(other_env, km1502_env)
other_env <- rbind(other_env, tn271_env)
other_env <- other_env[,c("time", "lat", "lon", "par", "salinity", "temp", "cruise")]
other_env <- rbind(other_env, km1923_env)

# KM1712, KM1713, KM1923, KM1502, SR1917, TN271, Gradients 4 nutrient data
# From Mercator-Pisces Biogeochemistry Weekly and Daily Forecast
km1712_nutr <- read.csv("~/Downloads/paper/raw_data/KM1712_nutr.csv")
km1712_nutr$time <- as.POSIXct(km1712_nutr$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
km1712_nutr$cruise <- "KM1712"

km1713_nutr <- read.csv("~/Downloads/paper/raw_data/KM1713_nutr.csv")
km1713_nutr$time <- as.POSIXct(km1713_nutr$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
km1713_nutr$cruise <- "KM1713"
km1713_nutr$day <- NULL
km1713_nutr$X <- NULL

km1502_nutr <- read.csv("~/Downloads/paper/raw_data/KM1502_nutr.csv")
km1502_nutr$time <- as.POSIXct(km1502_nutr$time, format = "%Y-%m-%dT%H:%M", tz = "GMT")
km1502_nutr$cruise <- "KM1502"

other_nutr <- rbind(km1712_nutr, km1713_nutr)
other_nutr <- rbind(other_nutr, km1502_nutr)

other_nutr <- subset(other_nutr, Fe < 2)
other_nutr <- subset(other_nutr, PP < 10)

nutr_lat <- group_by(other_nutr, lat, lon, cruise)
other_nutr <- dplyr::summarize_all(nutr_lat, mean)
other_nutr <- as.data.frame(other_nutr)
other_nutr$time <- as.POSIXct(other_nutr$time, format = "%Y-%m-%d", tz = "GMT")

other_nutr <- other_nutr[,c('lat','lon', 'time', 'NO3', 'PO4', 'Si', 'Fe', 'PP', 'CHL')]

############ combine data into one dataset

# get hourly averages for TN398 data
tn398_sfl$hour <- substr(tn398_sfl$time, 1, 13)
tn398_sfl$hour <- as.POSIXct(tn398_sfl$hour, format = "%Y-%m-%d %H", tz = "GMT")
sfl_hourly <- group_by(tn398_sfl, hour, pop)
sfl_hourly <- dplyr::summarize_all(sfl_hourly, mean)
sfl_hourly <- as.data.frame(sfl_hourly)
sfl_hourly$hour <- as.POSIXct(sfl_hourly$hour, format = "%Y-%m-%d %H", tz = "GMT")

tn398_env$hour <- substr(tn398_env$time, 1, 13)
tn398_env$hour <- as.POSIXct(tn398_env$hour, format = "%Y-%m-%d %H", tz = "GMT")
tn398_env_hr <- group_by(tn398_env, hour)
tn398_env_hr <- dplyr::summarize_all(tn398_env_hr, mean)
tn398_env_hr <- as.data.frame(tn398_env_hr)
tn398_env_hr$hour <- as.POSIXct(tn398_env_hr$hour, format = "%Y-%m-%d %H", tz = "GMT")

# combine environmental and phyto TN398 data
tn398 <- plyr::join(sfl_hourly, tn398_env_hr, by = "hour")
tn398$cruise <- "TN398"
tn398$LON <- round(tn398$LON, digits = 2)
tn398_nutr$Longitude <- round(tn398_nutr$Longitude, digits = 2)
tn398_nutr$LON <- tn398_nutr$Longitude
tn398 <- merge(tn398, tn398_nutr, by = "LON", all = TRUE)
tn398$NO3_NO2 <- tn398$NO3 + tn398$NO2

# get hourly averages for Gradients 1, 2, 3 data
grad_sfl$hour <- substr(grad_sfl$time, 1, 13)
grad_sfl$hour <- as.POSIXct(grad_sfl$hour, format = "%Y-%m-%d %H", tz = "GMT")
sfl_hourly <- group_by(grad_sfl, hour, pop, cruise)
sfl_hourly <- dplyr::summarize_all(sfl_hourly, mean)
sfl_hourly <- as.data.frame(sfl_hourly)
sfl_hourly$hour <- as.POSIXct(sfl_hourly$hour, format = "%Y-%m-%d %H", tz = "GMT")

grad_env$hour <- substr(grad_env$time, 1, 13)
grad_env$hour <- as.POSIXct(grad_env$hour, format = "%Y-%m-%d %H", tz = "GMT")
env_hourly <- group_by(grad_env, hour, cruise)
env_hourly <- dplyr::summarize_all(env_hourly, mean)
env_hourly <- as.data.frame(env_hourly)
env_hourly$hour <- as.POSIXct(env_hourly$hour, format = "%Y-%m-%d %H", tz = "GMT")

hourly <- merge(sfl_hourly, env_hourly, by = "hour")
hourly <- hourly[, c("hour", "pop", "cruise.x", "lat.x", "lon.x", "abundance", "diam", "Qc", "biomass", "par", "salinity", "temp")]
colnames(hourly) <- c("hour", "pop", "cruise", "lat", "lon", "abundance", "diam", "Qc", "biomass", "par", "salinity", "temp")

grad_nutr$hour <- substr(grad_nutr$time, 1, 13)
grad_nutr$hour <- as.POSIXct(grad_nutr$hour, format = "%Y-%m-%d %H", tz = "GMT")
grad_nutr_hr <- group_by(grad_nutr, hour)
grad_nutr_hr <- dplyr::summarize_all(grad_nutr_hr, mean)
grad_nutr_hr <- as.data.frame(grad_nutr_hr)
grad_nutr_hr$hour <- as.POSIXct(grad_nutr_hr$hour, format = "%Y-%m-%d %H", tz = "GMT")

# combine environmental and phyto Gradients data
grad <- plyr::join(hourly, grad_nutr_hr, by = "hour")

# get hourly averages for other data
other_sfl$hour <- substr(other_sfl$time, 1, 13)
other_sfl$hour <- as.POSIXct(other_sfl $hour, format = "%Y-%m-%d %H", tz = "GMT")
sfl_hourly <- group_by(other_sfl, hour, pop, cruise)
sfl_hourly <- dplyr::summarize_all(sfl_hourly, mean)
sfl_hourly <- as.data.frame(sfl_hourly)
sfl_hourly$hour <- as.POSIXct(sfl_hourly$hour, format = "%Y-%m-%d %H", tz = "GMT")

other_env$hour <- substr(other_env$time, 1, 13)
other_env$hour <- as.POSIXct(other_env$hour, format = "%Y-%m-%d %H", tz = "GMT")
env_hourly <- group_by(other_env, hour, cruise)
env_hourly <- dplyr::summarize_all(env_hourly, mean)
env_hourly <- as.data.frame(env_hourly)
env_hourly$hour <- as.POSIXct(env_hourly$hour, format = "%Y-%m-%d %H", tz = "GMT")

sfl_hourly$lat <- round(sfl_hourly$lat, digits = 1)
sfl_hourly$lon <- round(sfl_hourly$lon, digits = 1)
other <- plyr::join(sfl_hourly, other_nutr, by = c("lat", "lon"))

hourly <- merge(other, env_hourly, by = "hour")
hourly <- hourly[, c("hour", "pop", "cruise.x", "lat.x", "lon.x", "abundance", "diam", "Qc", "biomass", "par", "salinity", "temp", "NO3", "PO4", "Si", "Fe", "PP", "CHL")]
colnames(hourly) <- c("hour", "pop", "cruise", "lat", "lon", "abundance", "diam", "Qc", "biomass", "par", "salinity", "temp", "NO3_NO2", "PO4", "SiO4", "Fe", "PP", "chl")

# combining all data together
combo <- grad[,c("hour", "lat", "lon", "cruise", "pop", "abundance", "diam", "Qc", "biomass", "SiO4", "NO3_NO2", "PO4", "salinity", "temp", "par")]
to_add <- tn398[,c("hour", "lat", "lon", "cruise", "pop", "abundance", "diam", "Qc", "biomass", "Si.OH.4_uM", "NO3_NO2", "PO4", "SALINITY", "OCEAN TEMP", "PAR")]
colnames(to_add) <- c("hour", "lat", "lon", "cruise", "pop", "abundance", "diam", "Qc", "biomass", "SiO4", "NO3_NO2", "PO4", "salinity",  "temp", "par")
combo <- rbind(combo, to_add)
combo$Fe <- NA ##################### CHANGE THIS WHEN YOU GET GRADIENTS IRON DATA ######################
combo$PP <- NA
combo$chl <- NA
combo <- rbind(combo, hourly)
combo <- subset(combo, pop == "picoeuk" | pop == "prochloro" | pop == "synecho")

# exporting as csv
write.csv(combo, "~/Downloads/paper/all_data_PROGRESS.csv")

########

# test$time <- as.POSIXct(test$time, format = "%Y-%m-%d", tz = "GMT")
# test$day <- day(test$time)
# test <- subset(test, flag == 0)

# day <- unique(test$day)

# for (day_time in day){

# fcs <- subset(test, day == day_time)


# png(paste0("~/Downloads/paper","/test/", day_time, ".png"),width=8, height=6, unit="in", res=200)
# p <- list()
      # p[[1]] <- plot_vct_cytogram(fcs, "fsc_med","chl_med") + xlab("scatter") + ylab("red")
      # p[[2]] <- plot_cytogram(fcs, "fsc_med","chl_med") + xlab("scatter") + ylab("red")
      # p[[3]] <- plot_vct_cytogram(fcs, "fsc_med","pe_med") + xlab("scatter") + ylab("orange")
      # p[[4]] <- plot_cytogram(fcs, "fsc_med","pe_med") + xlab("scatter") + ylab("orange")
      # t <- textGrob(day_time, gp=gpar(fontsize = 20))
# grid.arrange(grobs = p, nrow = 2, top = t)
# dev.off()


# }


####### 
# modeled <- read.csv("~/Downloads/paper/raw_data/grad1_nutr_check.csv")
# modeled$time <- as.POSIXct(modeled$time, format = "%Y-%m-%d", tz = "GMT")

# test <- NULL
# test$time <- modeled$time
# test$NO3 <- modeled$NO3
# test <- as.data.frame(test)
# test$type <- "modeled"
# test$lat <- modeled$lat
# test <- subset(test, time > "2016-04-20" & time < "2016-05-04")

# add <- NULL
# add$time <- data_day$day
# add$NO3 <- data_day$NO3_NO2_mean
# add$type <- "observed"
# add <- as.data.frame(add)
# add$lat <- data_day$lat_mean
# add <- subset(add, time > "2016-04-20" & time < "2016-05-04")

# test <- rbind(test, add)

# png(paste0(folder, "nutr_test", ".png"),width=14, height=10, unit="in", res=200)
# p <- ggplot(data = test) +
	 # geom_point(aes(time, NO3, color = lat), alpha = 0.9, size=3, show.legend = T) +
	 # theme_bw() +
	 # theme(text = element_text(size = 20))
# print(p)
# dev.off()
