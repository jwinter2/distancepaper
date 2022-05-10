##################################################
Carbon Fixation Code
Written by Annette Hynes, adapted by Jordan Winter
Last edited 5/10/22
##################################################

library(popcycle)
library(tidyverse)
library(grid)
library(gridExtra)
library(suncalc)
library(googlesheets4)
library(viridis)
library(latex2exp)
library(openxlsx)
library(ggpubr)
library(lubridate)
library(marmap)
library(oce)
library(knitr)
library(broom)
library(ggh4x)
library(grDevices)

#################
### LOAD DATA ###
#################

save_path <- '~/Downloads/paper/'

file_name <- paste0(save_path, 'all_data_PROGRESS.csv')
all_data <- read.csv(file = file_name)

phyto_list <- c("prochloro", "synecho", "picoeuk", "croco")
all_data_phyto <- subset(all_data, pop %in% phyto_list)
all_data_phyto$pop <- factor(all_data_phyto$pop, levels = phyto_list)
rm(all_data)

t_data <- as.POSIXct(all_data_phyto$hour, format = '%Y-%m-%d %H:%M:%S', tz = 'GMT')
all_data_phyto$time <- t_data

all_data_phyto$year <- as.factor(lubridate::year(all_data_phyto$time))
all_data_phyto$month <- as.factor(lubridate::month(all_data_phyto$time))
all_data_phyto$day <- as.factor(lubridate::yday(all_data_phyto$time))
all_data_phyto$hour <- as.factor(lubridate::hour(all_data_phyto$time))
all_data_phyto$date <- lubridate::as_date(all_data_phyto$time)

group.colors <- c(prochloro = viridis::viridis(4)[1],
        synecho = viridis::viridis(4)[2], picoeuk = viridis::viridis(4)[3])

col_seaflow <- viridis::viridis(9)[5]
col_HOT <- '#FFFFFFFF'

prog.colors <- c(SeaFlow = '#21908CAA', HOT = '#FFFFFFFF')

####################
### HOURLY MEANS ###
####################

data_in <- all_data_phyto[, c("lat", "lon", "cruise", "pop", 'abundance', 'diam', 'Qc', 'biomass', 'date', 'time')]

# Mid-month sunrise and sunset for plotting
midmonth <- as.Date(paste0(lubridate::year(data_in$time), '-', lubridate::month(data_in$time), '-15'), tz = 'GMT')

# get sunrise and sunset times using suncalc package

sun_mid <- NULL
sun_mid$date <- NULL
sun_mid$sunrise <- NULL
sun_mid$sunset <- NULL
for (i in 1:nrow(data_in)){
sun_mid_row <- suncalc::getSunlightTimes(date = data_in$date[i], lat = data_in$lat[i], lon = data_in$lon[i], keep = c('sunrise', 'sunset'), tz = "GMT")
sun_mid$date[i] <- as.character(sun_mid_row$date)
sun_mid$sunrise[i] <- as.character(sun_mid_row$sunrise)
sun_mid$sunset[i] <- as.character(sun_mid_row$sunset)
}
sun_mid <- as.data.frame(sun_mid)
sun_mid$date <- as.POSIXct(sun_mid$date, format = '%Y-%m-%d', tz = 'GMT')
sun_mid$sunrise <- as.POSIXct(sun_mid$sunrise, format = '%Y-%m-%d %H:%M:%S', tz = 'GMT')
sun_mid$sunset <- as.POSIXct(sun_mid$sunset, format = '%Y-%m-%d %H:%M:%S', tz = 'GMT')
sunrise_mid <- lubridate::hour(sun_mid$sunrise) + lubridate::minute(sun_mid$sunrise)/60

data_in$time_seconds

sunrise_seconds <- sunrise_mid*3600
suntime <- data_in$time - sunrise_seconds
data_in$sundate <- lubridate::as_date(suntime)
data_in$sunhour <- lubridate::hour(suntime)
data_in$sundate <- as.factor(data_in$sundate)

data_in$pop <- factor(data_in$pop, levels = c('prochloro', 'synecho', 'picoeuk'))
data_in <- subset(data_in, is.na(sundate) == F)

# Get grouped daily means to detrend data and get fold change

all_datasun <- data_in %>%
    group_by(pop, sundate) %>%
    mutate(abundance_detrend = abundance - mean(abundance, na.rm = TRUE), diam_detrend = diam - mean(diam, na.rm = TRUE), Qc_detrend = Qc - mean(Qc, na.rm = TRUE),
        abundance_fold = if(dplyr::first(sunhour == 0))             # Use the dawn point.  If it doesn't exist, use the last hour (just before the next dawn)
            (abundance - dplyr::first(abundance))/first(abundance)
            else (abundance - dplyr::last(abundance))/last(abundance),
        diam_fold = if(dplyr::first(sunhour == 0))
            (diam - dplyr::first(diam))/first(diam)
            else (diam - dplyr::last(diam))/last(diam),
        Qc_fold = if(dplyr::first(sunhour == 0))
            (Qc - dplyr::first(Qc))/first(Qc)
            else (Qc - dplyr::last(Qc))/last(Qc))


#########################
### C FIXATION: DAILY ###
#########################

# Take only hourly means from daylight hours.  Use data table "robust" with Qc outlier cruises removed and abundance >= 0.02.
# Remove days with less than 6 hours between first and last daily point.
all_data_hr_robust <- data_in
sun <- NULL
sun$date <- NULL
sun$sunrise <- NULL
sun$sunset <- NULL
for (i in 1:nrow(data_in)){
sun_mid_row <- suncalc::getSunlightTimes(date = all_data_hr_robust$date[i], lat = all_data_hr_robust$lat[i], lon = all_data_hr_robust$lon[i], keep = c('sunrise', 'sunset'), tz = "GMT")
sun$date[i] <- as.character(sun_mid_row$date)
sun$sunrise[i] <- as.character(sun_mid_row$sunrise)
sun$sunset[i] <- as.character(sun_mid_row$sunset)
}
sun <- as.data.frame(sun)
all_data_hr_robust$date <- as.character(all_data_hr_robust$date)
sun$date <- as.character(sun$date)
all_data_hr_robust <- merge(sun, all_data_hr_robust, by = "date")

all_data_hr_robust$sunrise <- as.POSIXct(all_data_hr_robust$sunrise, format = '%Y-%m-%d %H:%M:%S', tz = 'GMT')
all_data_hr_robust$sunset <- as.POSIXct(all_data_hr_robust$sunset, format = '%Y-%m-%d %H:%M:%S', tz = 'GMT')

# Remove low abundance points first
all_data_hr_robust <- subset(all_data_hr_robust, abundance > 0.02) # 0.048 = About 30 cells per 3-min file

all_data_hr_robust$light <- "night"
ind_day <- which(all_data_hr_robust$time < all_data_hr_robust$sunset & all_data_hr_robust$time > all_data_hr_robust$sunrise)
all_data_hr_robust$light[ind_day] <- "day"
all_data_hr_robust$daylength <- as.numeric(difftime(all_data_hr_robust$sunset, all_data_hr_robust$sunrise, units = "hours"))
all_data_light <- subset(as.data.frame(all_data_hr_robust), light == "day")
all_data_light$hour <- as.numeric(lubridate::hour(all_data_light$time))    # Otherwise the linear model treats it as categorical

# Linear regression to estimate rate of C fixation
data_lm_day <- all_data_light %>%
    dplyr::group_by(date, pop, cruise) %>%
    tidyr::nest() %>%
    dplyr::mutate(
        model = purrr::map(data, ~lm(Qc ~ hour, data = .)),
        tidied = purrr::map(model, broom::tidy)
    ) %>%
    tidyr::unnest(tidied)

data_lm_day_df <- as.data.frame(data_lm_day[data_lm_day$term == "hour", c("pop", 'date', 'estimate', 'std.error', 'cruise', 'p.value')])

Qc_mean <- all_data_light %>%
    dplyr::group_by(date, pop, cruise) %>%
    dplyr::summarize(Qc_mean = mean(Qc, na.rm = TRUE), data_day = last(hour) - first(hour))

data_lm_day_df <- merge(Qc_mean, data_lm_day_df)
data_lm_day_df$date <- as.POSIXct(data_lm_day_df$date, format = '%Y-%m-%d', tz = 'GMT')

# Quality control
data_lm_day_df$p.value[which(data_lm_day_df$p.value >= 0.05)] <- NA # Remove insignificant linear regressions
data_lm_day_df <- data_lm_day_df[!is.na(data_lm_day_df$p.value), ] # Remove flagged p.values

data_lm_day_df$estimate[which(data_lm_day_df$estimate < 0)] <- NA   # Remove negative C fixation estimates
data_lm_day_df <- data_lm_day_df[!is.na(data_lm_day_df$estimate), ]  # Remove lines with NA C fixation estimate
data_lm_day_df$data_day[which(data_lm_day_df$data_day < 6)] <- NA    # Tag days when first and last data points are less than 6 hours apart
data_lm_day_df <- data_lm_day_df[!is.na(data_lm_day_df$data_day), ]  # Remove lines with short data days

data_lm_day_df$Cfix_norm <- data_lm_day_df$estimate/data_lm_day_df$Qc_mean
data_lm_day_df$Cfix_se <- data_lm_day_df$std.error/data_lm_day_df$Qc_mean
data_lm_day_df$Cmax <- data_lm_day_df$Cfix_norm + data_lm_day_df$Cfix_se
data_lm_day_df$Cmin <- data_lm_day_df$Cfix_norm - data_lm_day_df$Cfix_se

data_lm_day_df <- subset(data_lm_day_df, is.na(pop) == F)

# make figures

fig_name <- paste0(save_path, "Qc_fixation_norm_date.pdf")
pdf(fig_name, width = 16, height = 10)
g <- ggplot(data_lm_day_df) +
    #ggplot2::geom_linerange(ggplot2::aes(x = date, ymin = Cmin, ymax = Cmax), color = "grey") +
    geom_point(aes_string(x = "date", y = "Cfix_norm", fill = "pop"), pch = 21, size = 3, alpha = 0.5) +
    scale_fill_manual(values = group.colors) +
    theme_bw(base_size = 10) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = "Date", y = latex2exp::TeX('Net C fixation (fraction cellular carbon hr$^{-1}$)'), title = "Daily Average Normalized Carbon Fixation") +
    facet_wrap(pop ~ cruise, nrow = 3, scales = 'free')
print(g)
dev.off()

fig_name <- paste0(save_path, "Qc_fixation_date.pdf")
pdf(fig_name, width = 18, height = 10)
g <- ggplot(data_lm_day_df) +
    #ggplot2::geom_linerange(ggplot2::aes(x = date, ymin = Cmin, ymax = Cmax), color = "grey") +
    geom_point(aes_string(x = "date", y = "estimate", fill = "pop"), pch = 21, size = 3, alpha = 0.5) +
    scale_fill_manual(values = group.colors) +
    theme_bw(base_size = 10) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = "Date", y = latex2exp::TeX('Net C fixation (pg C cell$^{-1}$ hr$^{-1}$)'), title = "Daily Average Carbon Fixation Rates") +
    facet_wrap(pop ~ cruise, nrow = 3, scales = 'free')
print(g)
dev.off()

fig_name <- paste0(save_path, "Qc_mean_date.pdf")
pdf(fig_name, width = 16, height = 10)
g <- ggplot(data_lm_day_df) +
    #ggplot2::geom_linerange(ggplot2::aes(x = date, ymin = Cmin, ymax = Cmax), color = "grey") +
    geom_point(aes_string(x = "date", y = "Qc_mean", fill = "pop"), pch = 21, size = 3, alpha = 0.5) +
    scale_fill_manual(values = group.colors) +
    theme_bw(base_size = 10) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = "Date", y = latex2exp::TeX('Mean Qc (pg C cell$^{-1}$)'), title = "Daily Average Carbon Quota") +
    facet_wrap(pop ~ cruise, nrow = 3, scales = 'free')
print(g)
dev.off()

data_lm_day_df$month <- lubridate::month(data_lm_day_df$date)
data_lm_day_df$day <- lubridate::day(data_lm_day_df$date)
data_lm_day_df$doy <- lubridate::yday(data_lm_day_df$date)
data_lm_day_df$year <- factor(lubridate::year(data_lm_day_df$date))

fig_name <- paste0("~/Downloads/paper/other_figures/", "Cspec_growthrate.pdf")
pdf(fig_name, width = 10, height = 8)
g2 <- ggplot(data_lm_day_df, aes_string(x = 'Cfix_norm')) +
        geom_histogram(alpha = 0.5, fill = "grey30") +
        theme_bw(base_size = 10) +
        facet_wrap(vars(pop), ncol = 1) +
        ggplot2::labs(y = 'No. days', x = "Carbon Specific Growth Rate")
print(g2)
dev.off()

# export csv file

write.csv(data_lm_day_df, "~/Downloads/paper/Cfix.csv")

#### combine Cfix data with temp and sal for plotting

combo_no_croco <- read.csv("~/Downloads/paper/all_data_PROGRESS.csv")
combo_no_croco <- subset(data, pop != "croco")
combo_no_croco$hour <- as.POSIXct(combo_no_croco$hour, format = "%Y-%m-%d %H", tz = "GMT")
combo_no_croco$date <- substr(combo_no_croco$hour, 1, 10)
combo_no_croco$date <- as.character(combo_no_croco$date)

to_add <- combo_no_croco[, c("date", "salinity", "temp", "SiO4", "NO3_NO2", "PO4")]
grouping <- group_by(to_add, date)
to_add <- dplyr::summarize_all(grouping, mean, na.rm = T)
to_add <- as.data.frame(to_add)

data_lm_day_df$date <- as.character(data_lm_day_df$date)

total_qc <- merge(data_lm_day_df, to_add, by = "date")
total_qc <- subset(total_qc, temp < 40)

png(paste0(folder, "cfix_temp_cruise.png"),width=18, height=10, unit="in", res=200)
g <- ggplot() +
	 geom_point(data = total_qc, aes(temp, Cfix_norm, color = cruise), size=3, show.legend = T) +
	 labs(x='Temperature (C)', y= 'Carbon Specific Growth Rate') +
	 theme_bw() +
	 theme(text = element_text(size = 15)) +
	 scale_color_manual(values = viridis(10)) +
	 ggtitle("Carbon Fixation Rates and Temperature") +
	 facet_wrap(pop ~ ., nrow = 1, scales = "free_y")
print(g)
dev.off()

png(paste0(folder, "cfix_sal_cruise.png"),width=18, height=10, unit="in", res=200)
g <- ggplot() +
	 geom_point(data = total_qc, aes(salinity, Cfix_norm, color = cruise), size=3, show.legend = T) +
	 labs(x='Salinity (PSU)', y= 'Carbon Specific Growth Rate') +
	 theme_bw() +
	 theme(text = element_text(size = 15)) +
	 scale_color_manual(values = viridis(10)) +
	 ggtitle("Carbon Fixation Rates and Salinity") +
	 facet_wrap(pop ~ ., nrow = 1, scales = "free_y")
print(g)
dev.off()

png(paste0(folder, "cfix_no3_cruise.png"),width=18, height=10, unit="in", res=200)
g <- ggplot() +
	 geom_point(data = total_qc, aes(NO3_NO2, Cfix_norm, color = cruise), size=3, show.legend = T) +
	 labs(x='Nitrate and Nitrite Concentration (ug/L)', y= 'Carbon Specific Growth Rate') +
	 theme_bw() +
	 theme(text = element_text(size = 15)) +
	 scale_color_manual(values = viridis(10)) +
	 ggtitle("Carbon Fixation Rates and Nitrate") +
	 facet_wrap(pop ~ ., nrow = 1, scales = "free_y")
print(g)
dev.off()

png(paste0(folder, "cfix_po4_cruise.png"),width=18, height=10, unit="in", res=200)
g <- ggplot() +
	 geom_point(data = total_qc, aes(PO4, Cfix_norm, color = cruise), size=3, show.legend = T) +
	 labs(x='Phosphate Concentration (ug/L)', y= 'Carbon Specific Growth Rate') +
	 theme_bw() +
	 theme(text = element_text(size = 15)) +
	 scale_color_manual(values = viridis(10)) +
	 ggtitle("Carbon Fixation Rates and Phosphate") +
	 facet_wrap(pop ~ ., nrow = 1, scales = "free_y")
print(g)
dev.off()

##################
### Check CFIX ###
##################

# Each day, all phytoplankton together, to see if the rates appear to be reasonable
data_lm_day_df$pop <- factor(data_lm_day_df$pop, levels = names(group.colors))  # match population order

day_list <- unique(data_lm_day_df$date)  # The days for which we have C fixation rates
for (i in seq(1, length(day_list))){
    siku <- day_list[i]
    print(siku)
    leo <- subset(all_data_hr_robust, sundate == as.character(siku))
    PP_leo <- subset(data_lm_day_df, date == (siku))

    g1 <- ggplot(leo, aes_string(x = 'sunhour', y = 'Qc')) +
        geom_rect(aes(xmin = daylength[1], xmax = 24, ymin = -Inf, ymax = Inf), fill = "grey", color = 'grey') +
        geom_line(color = prog.colors['SeaFlow'], size = 2) +
        geom_abline(data = PP_leo, aes(intercept = min(Qc_mean), slope = estimate), linetype = "dashed", size = 1) +
        theme_bw(base_size = 20) +
        facet_wrap(vars(pop), ncol = 1, scales = 'free_y', drop = FALSE) +
        ggplot2::labs(y = unname(latex2exp::TeX('Carbon Quota (pg C cell$^{-1}$)')), x = 'Hours since Dawn')

    g2 <- ggplot(data_lm_day_df, aes_string(x = 'estimate')) +
        geom_histogram(alpha = 0.5, fill = "grey30") +
        geom_vline(data = PP_leo, aes(xintercept = estimate), size = 1) +
        theme_bw(base_size = 20) +
        ggplot2::theme(legend.position = 'none') +
        facet_wrap(vars(pop), ncol = 1, scales = 'free_x') +
        ggplot2::labs(y = 'No. days', x = unname(latex2exp::TeX('C fix (pg C cell$^{-1}$ hr$^{-1}$)')))

    g3 <- ggplot(data_lm_day_df, aes_string(x = 'Qc_mean')) +
        geom_histogram(alpha = 0.5, fill = "grey30") +
        geom_vline(data = PP_leo, aes(xintercept = Qc_mean), size = 1) +
        theme_bw(base_size = 20) +
        ggplot2::theme(legend.position = 'none') +
        facet_wrap(vars(pop), ncol = 1, scales = 'free_x') +
        ggplot2::labs(y = 'No. days', x = unname(latex2exp::TeX('Daytime Mean Qc (pg C cell$^{-1}$)')))

    g4 <- ggplot(data_lm_day_df, aes_string(x = 'Cfix_norm')) +
        geom_histogram(alpha = 0.5, fill = "grey30") +
        geom_vline(data = PP_leo, aes(xintercept = Cfix_norm), size = 1) +
        theme_bw(base_size = 20) +
        ggplot2::theme(legend.position = 'none') +
        facet_wrap(vars(pop), ncol = 1) +
        ggplot2::labs(y = 'No. days', x = unname(latex2exp::TeX('C fix norm (hr$^{-1}$)')))

    t <- textGrob(siku, gp = gpar(fontsize=30))

    fig_file <- paste0(save_path, 'Cfix/Cfix_by_date_', siku, '.png')
    png(fig_file, width = 1500, height = 1200)
    grid.arrange(g1, g2, g3, g4, ncol = 4, top = t)
    dev.off()

}   # end day loop
