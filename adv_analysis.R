# This code analyzes the data from the Nortek Vector ADV.
# taken from tke_reader.r, currently working on this code.

# to install the needed packages, run the following lines:
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("lubridate")
library(dplyr)
library(ggplot2)
library(lubridate)

setwd("c:/Users/duquesne/Documents/nortek/data") # lab laptop
setwd("/Users/davidkahler/Documents/Hydrology_and_WRM/river_and_lake_mixing/ADV_data/") # David's computer
fh <- "528AR214" # filename header
fn_sen <- paste(fh, "sen", sep = ".")
sen <- read.table(fn_sen, header = FALSE, sep = "", dec = ".")
sen <- sen %>% rename(mon = V1, day = V2, yea = V3, hou = V4, mnt = V5, sec = V6, err = V7, sta = V8, bat = V9, ssp = V10, hed = V11, pit = V12, rol = V13, tmp = V14, a1 = V15, checksum = V16)
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
d <- 24*3600*as.numeric(as.Date(paste(sen$yea[1], sen$mon[1], sen$day[1], sep="-"), origin="1970-01-01")) # number of seconds that gives the day
h <- sen$sec[1]+60*sen$mnt[1]+3600*sen$hou[1] # time in seconds
starttime <- as_datetime(d + h) # lubridate datetime for the start of the data

fn_dat <- paste(fh, "dat", sep = ".")
dat <- read.table(fn_dat, header = FALSE, sep = "", dec = ".")
# or
# x <- file.choose()
# dat <- read.table(x, header = FALSE, sep = "", dec = ".")
dat <- dat %>% rename(burst = V1, ensemble = V2, u = V3, v = V4, w = V5, amp1 = V6, amp2 = V7, amp3 = V8, snr1 = V9, snr2 = V10, snr3 = V11, corr1 = V12, corr2 = V13, corr3 = V14, p_dbar = V15, a1 = V16, a2 = V17, checksum = V18)
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
#pck = read.table("MON103.pck", header = FALSE, sep = "", dec = ".")
#vhd = read.table("MON103.vhd", header = FALSE, sep = "", dec = ".")
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
plot(dat$time,(1e5*dat$p_dbar), type = "l",ylab = "Pressure (Pa)", xlab = "Time (s)")
lines(c(min(dat$time),max(dat$time)),c(101325,101325)) # Places a line at what should be the surface of the water

## Figure out an estimate of pressure:
# temperature <- array(NA, dim = nrow(dat))
# for (i in 1:nrow(dat)) {
#       if (dat$p_dbar>45000) {
#             temperature[i] <- # WE NEED TO PULL TEMP DATA FROM SEN...
#       }
# }
atmos <- mean(1e5*dat$p_dbar[1:10]) # Pa, to subtract atmospheric pressure
dat$depth <- -(1e5*dat$p_dbar - atmos)/(9.81*997)
plot(hms::as_hms(dat$time),(dat$depth), type = "l",ylab = "Depth (m)", xlab = "Time (s, from midnight)")

par(mfrow = c(3,1), mar = c(4,4,1,1))
plot(hms::as_hms(dat$time),dat$u, ylim = c(-1, 1), type = "l",ylab = "u (m/s)", xlab = "")
plot(hms::as_hms(dat$time),dat$v, ylim = c(-1, 1), type = "l",ylab = "v (m/s)", xlab = "")
plot(hms::as_hms(dat$time),dat$w, ylim = c(-1, 1), type = "l",ylab = "w (m/s)", xlab = "Time (s, from midnight)")

# ANALYSIS BY AVERAGING WINDOW:
u <- array(NA, dim = c(ceiling(nrow(dat)/bar),bar))
v <- u
w <- v
time <- array(-9999, dim = c(ceiling(nrow(dat)/bar))) # will identify the start of every averaged window.
u_ave <- array(0, dim = c(nrow(u),2))
v_ave <- u_ave
w_ave <- u_ave
u_prime <- u
v_prime <- v
w_prime <- w
uu <- array(0, dim = c(nrow(u)))
vv <- uu
ww <- uu
uv <- uu
uw <- uu
vw <- uu
for (i in 1:(nrow(u))) {
  time[i] <- dat$time[(bar*(i-1)) + 1] # this is the start time of the averaging window at time step, i
  for (j in 1:bar) { # this will cycle over each averaging window
    dat_index <- (bar*(i-1)) + j
    if (is.na(dat$u[dat_index])==FALSE) {
      if (dat$checksum[dat_index]>0) {
        print(paste0("error at ", dat$ensemble[dat_index]))
      }
      u[i,j] <- dat$u[dat_index] # organize data
      u_ave[i,2] = u_ave[i,2] + 1
      v[i,j] <- dat$v[dat_index]
      v_ave[i,2] = v_ave[i,2] + 1
      w[i,j] <- dat$w[dat_index]
      w_ave[i,2] = w_ave[i,2] + 1
    }
  }
  u_ave[i,1] <- mean(u[i,], na.rm = TRUE)
  v_ave[i,1] <- mean(v[i,], na.rm = TRUE)
  w_ave[i,1] <- mean(w[i,], na.rm = TRUE)
  for (j in 1:bar) {
    u_prime[i,j] <- u[i,j] - u_ave[i,1]
    v_prime[i,j] <- v[i,j] - v_ave[i,1]
    w_prime[i,j] <- w[i,j] - w_ave[i,1]
  }
  uu[i] <- mean((u_prime[i,]^2))
  vv[i] <- mean((v_prime[i,]^2))
  ww[i] <- mean((w_prime[i,]^2))
  uv[i] <- mean((u_prime[i,]*v_prime[i,]))
  uw[i] <- mean((u_prime[i,]*w_prime[i,]))
  vw[i] <- mean((v_prime[i,]*w_prime[i,]))
}

par(mfrow = c(3,1), mar = c(4,4,1,1))
# REM: x-axis is datetime range
plot(hms::as_hms(time),u_ave[,1], ylim = c(-1, 1), type = "l",ylab = "u (m/s)", xlab = "")
plot(hms::as_hms(time),v_ave[,1], ylim = c(-1, 1), type = "l",ylab = "v (m/s)", xlab = "")
plot(hms::as_hms(time),w_ave[,1], ylim = c(-1, 1), type = "l",ylab = "w (m/s)", xlab = "Time (s)")

# to zoom in on an area of interest, find indices:
start <- as.numeric(ymd_hms("2021-05-28 16:57:00")) # Enter start time here as "YYYY-MM-DD HH:MM:SS" in 24-hour time
end <- as.numeric(ymd_hms("2021-05-28 16:59:00")) # End time, same format
diff_s <- abs(time-start) # finds the difference between the time entries and start time (in case we don't hit it exactly)
diff_e <- abs(time-end)
m_s <- 10*bar_s # allocate variable for minimum finding.  larger than what will be found
m_e <- m_s
s <- NA
e <- NA
for (i in 1:(length(time))) {
      if (diff_s[i] < m_s) {
            m_s <- diff_s[i]
            s <- i # index of start
      }
      if (diff_e[i] <= m_e) {
            m_e <- diff_e[i]
            e <- i # index of end
      }
}

par(mfrow = c(3,1), mar = c(4,4,1,1))
plot(hms::as_hms(time[s:e]), uu[s:e], ylim = c(0, 1), type = "l", ylab = "uu", xlab = "")
plot(hms::as_hms(time[s:e]), vv[s:e], ylim = c(0, 1), type = "l", ylab = "vv", xlab = "")
plot(hms::as_hms(time[s:e]), ww[s:e], ylim = c(0, 1), type = "l", ylab = "ww", xlab = "Time (s)")

par(mfrow = c(3,1), mar = c(4,4,1,1))
plot(uv, ylim = c(-0.3, 0.3), type = "l", ylab = "uv", xlab = "")
plot(uw, ylim = c(-0.3, 0.3), type = "l", ylab = "uw", xlab = "")
plot(vw, ylim = c(-0.3, 0.3), type = "l", ylab = "vw", xlab = "Time (s)")

# Spectra
par(mfrow = c(3,1), mar = c(4,4,2,2))
hist(uv[start:stop,1], breaks = c(-10000,-1.5,-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1.5,10000), xlim = c(-1.5,1.5), ylab = "u'v'", xlab = "", main = "")
hist(uw[start:stop,1], breaks = c(-10000,-1.5,-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1.5,10000), xlim = c(-1.5,1.5), ylab = "u'w'", xlab = "", main = "")
hist(vw[start:stop,1], breaks = c(-10000,-1.5,-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1.5,10000), xlim = c(-1.5,1.5), ylab = "v'w'", xlab = "", main = "")



