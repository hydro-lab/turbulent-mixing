# This code pulls data from the Nortek Vector ADV and preps the endpoints for the profiling stops.  The 
# code will take both the .sen file with time and .dat file with pressure to determine the depth and 
# times.

library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(latex2exp)

#setwd("c:/Users/duquesne/Documents/nortek/data") # lab laptop
setwd("/Users/davidkahler/Documents/Hydrology_and_WRM/river_and_lake_mixing/ADV_data/") # David's computer
fh <- "109MON18" # filename header
fn_sen <- paste(fh, "sen", sep = ".")
sen <- read_table(fn_sen, col_names = FALSE, col_types = "nnnnnnnnnnnnnnnn")
sen <- sen %>%
      rename(mon = X1, day = X2, yea = X3, hou = X4, mnt = X5, sec = X6, err = X7, sta = X8, bat = X9, ssp = X10, hed = X11, pit = X12, rol = X13, tmp = X14, a1 = X15, checksum = X16) %>%
      mutate(dt = ymd_hms(paste(yea,mon,day,hou,mnt,sec)))
# 1   Month                            (1-12)
# 2   Day                              (1-31)
# 3   Year
# 4   Hour                             (0-23)
# 5   Minute                           (0-59)
# 6   Second                           (0-59)
# 7   Error code
# 8   Status code
# 9   Battery voltage                  (V)
# 10   Soundspeed                       (m/s)
# 11   Heading                          (degrees)
# 12   Pitch                            (degrees)
# 13   Roll                             (degrees)
# 14   Temperature                      (degrees C)
# 15   Analog input
# 16   Checksum                         (1=failed)

fn_dat <- paste(fh, "dat", sep = ".")
dat <- read_table(fn_dat, col_names = FALSE, col_types = "nndddnnnnnnnnnnnnn")
dat <- dat %>%
      rename(burst = V1, ensemble = V2, u = X3, v = X4, w = X5, amp1 = X6, amp2 = X7, amp3 = X8, snr1 = X9, snr2 = X10, snr3 = X11, corr1 = X12, corr2 = X13, corr3 = X14, p_dbar = X15, a1 = X16, a2 = X17, checksum = X18)
# 1   Burst counter
# 2   Ensemble counter                 (1-65536)
# 3   Velocity (Beam1|X|East)          (m/s)
# 4   Velocity (Beam2|Y|North)         (m/s)
# 5   Velocity (Beam3|Z|Up)            (m/s)
# 6   Amplitude (Beam1)                (counts)
# 7   Amplitude (Beam2)                (counts)
# 8   Amplitude (Beam3)                (counts)
# 9   SNR (Beam1)                      (dB)
# 10   SNR (Beam2)                      (dB)
# 11   SNR (Beam3)                      (dB)
# 12   Correlation (Beam1)              (%)
# 13   Correlation (Beam2)              (%)
# 14   Correlation (Beam3)              (%)
# 15   Pressure                         (dbar)        p_dbar
# 16   Analog input 1
# 17   Analog input 2
# 18   Checksum                         (1=failed)
sampling_rate = 64 # Hz, verify sampling rate in .hdr file under User setup

# Data check
records <- nrow(dat)
seconds <- nrow(sen)
top_samples <- seconds*sampling_rate # this is the number of records that there would be if every second had all of the sampling rate's elements filled in, the value (records) should not exceed this.
numbers <- ((records<=top_samples)&&(records>(top_samples-(2*sampling_rate)))) # we allow twice the values to account for missing values at the first second and the last second.
error_checksum <- max(dat$checksum) # If there are any errors recorded, this will be = 1
par(mfrow = c(3,1), mar = c(4,4,1,1)) # mfrow=c(nrows, ncols), https://www.statmethods.net/advgraphs/layout.html
hist(dat$snr1, breaks = c(-100,0,5,10,15,20,25,30,35,40,45,50,55,60,100), xlim = c(0,60), ylab = "Beam 1", xlab = "", main = "")
hist(dat$snr2, breaks = c(-100,0,5,10,15,20,25,30,35,40,45,50,55,60,100), xlim = c(0,60), ylab = "Beam 2", xlab = "", main = "")
hist(dat$snr3, breaks = c(-100,0,5,10,15,20,25,30,35,40,45,50,55,60,100), xlim = c(0,60), ylab = "Beam 3", xlab = "Signal-to-Noise Ratio (dB)", main = "")

# bar is the averaging window for U_bar and to compute the deviations from the mean, easiest to express 
# as a multiple of the sampling rate, therefore measured in seconds: bar_s.  Can take non-integer values.
bar_s <- 15 # averaging window in seconds.  Should evaluate range from depth/mean representative velocity to entire period
bar <- round(bar_s * sampling_rate) # in indexed values [i], round() needed to ensure that it fits within the dataset
# gives decimal time for each data record
dat$time <- starttime + (c(0:(nrow(dat)-1)))/sampling_rate

## Examine data:
par(mfrow = c(1,1))
plot(dat$time,(1e4*dat$p_dbar), type = "l",ylab = "Pressure (Pa)", xlab = "Time (s)")
lines(c(min(dat$time),max(dat$time)),c(101325,101325)) # Places a line at what should be the surface of the water

## Figure out an estimate of pressure:
# temperature <- array(NA, dim = nrow(dat))
# for (i in 1:nrow(dat)) {
#       if (dat$p_dbar>45000) {
#             temperature[i] <- # WE NEED TO PULL TEMP DATA FROM SEN...
#       }
# }
atmos <- mean(1e4*dat$p_dbar[1:10]) # Pa, to subtract atmospheric pressure
dat$depth <- -(1e4*dat$p_dbar - atmos)/(9.81*997)
#xlim = c(ymd_hms("2021-06-29 15:09:00",ymd_hms("2021-06-29 15:11:00")))
par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(dat$time,dat$depth, type = "l", ylim = c(-6,0), ylab = "Depth (m)", xlab = "Time")

par(mfrow = c(3,1), mar = c(4,4,1,1))
plot(hms::as_hms(dat$time),dat$u, ylim = c(-10, 10), type = "l",ylab = "u (m/s)", xlab = "")
plot(hms::as_hms(dat$time),dat$v, ylim = c(-10, 10), type = "l",ylab = "v (m/s)", xlab = "")
plot(hms::as_hms(dat$time),dat$w, ylim = c(-10, 10), type = "l",ylab = "w (m/s)", xlab = "Time (s, from midnight)")
