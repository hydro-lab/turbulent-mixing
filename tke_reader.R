# This code analyzes the data from the Nortek Vector ADV.

# to install the needed packages, run the following lines:
#install.packages("dplyr")
#install.packages("ggplot2")
library(dplyr)
library(ggplot2)

setwd("/Users/davidkahler/Documents/Hydrology_and_WRM/river_and_lake_mixing/ADV_data/")
fh <- "chartiers1" # filename header
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
# 15   Pressure                         (dbar)
# 16   Analog input 1
# 17   Analog input 2
# 18   Checksum                         (1=failed)
#pck = read.table("MON103.pck", header = FALSE, sep = "", dec = ".")
#vhd = read.table("MON103.vhd", header = FALSE, sep = "", dec = ".")
sampling_rate = 64 # Hz, verify sampling rate in .hdr file under User setup

# Data check
records <- nrow(dat)
seconds <- nrow(sen)
top_samples <- seconds*sampling_rate # this is the number of records that there would be if every second had all of the sampling rate's elements filled in, the value should not exceed this.
numbers <- ((records<=top_samples)&&(records>(top_samples-(2*sampling_rate)))) # we allow twice the values to account for missing values at the first second and the last second.
error_checksum <- max(dat$checksum) # If there are any errors recorded, this will be = 1
par(mfrow = c(3,1), mar = c(6,6,3,3)) # mfrow=c(nrows, ncols), https://www.statmethods.net/advgraphs/layout.html
hist(dat$snr1, breaks = c(-100,0,5,10,15,20,25,30,35,40,45,50,55,60,100), xlim = c(0,60), ylab = "Beam 1", xlab = "", main = "")
hist(dat$snr2, breaks = c(-100,0,5,10,15,20,25,30,35,40,45,50,55,60,100), xlim = c(0,60), ylab = "Beam 2", xlab = "", main = "")
hist(dat$snr3, breaks = c(-100,0,5,10,15,20,25,30,35,40,45,50,55,60,100), xlim = c(0,60), ylab = "Beam 3", xlab = "Signal-to-Noise Ratio (dB)", main = "")

datetime <- array(-9999, dim = c(nrow(sen),2)) # col 1: days (day 1 is 01 Jan 2020), col 2: seconds of the day
u <- array(-9999, dim = c(nrow(sen),sampling_rate))
v <- array(-9999, dim = c(nrow(sen),sampling_rate))
w <- array(-9999, dim = c(nrow(sen),sampling_rate))
u_ave <- array(0, dim = c(nrow(sen),2))
v_ave <- array(0, dim = c(nrow(sen),2))
w_ave <- array(0, dim = c(nrow(sen),2))
u_prime <- u
v_prime <- v
w_prime <- w
uu <- array(-9999, dim = c(nrow(sen),1))
vv <- uu
ww <- uu
uv <- uu
uw <- uu
vw <- uu
for (i in 1:nrow(sen)) {
      datetime[i,1] <- ( as.numeric(as.Date(paste(sen$yea[i], sen$mon[i], sen$day[i], sep = " "), format = "%Y %m %d")) - as.numeric(as.Date("2019 12 31", format = "%Y %m %d")) )
      datetime[i,2] <- sen$sec[i] + (sen$mnt[i])*60 + (sen$hou[i])*3600
      for (j in 1:sampling_rate) {
            dat_index <- (sampling_rate*(i-1)) + j
            if (is.na(dat$u[dat_index])==FALSE) {
                  if (dat$checksum[dat_index]>0) {
                        print(paste0("error at ", dat$ensemble[dat_index]))
                  }
                  u[i,j] <- dat$u[dat_index]
                  u_ave[i,2] = u_ave[i,2] + 1
                  v[i,j] <- dat$v[dat_index]
                  v_ave[i,2] = v_ave[i,2] + 1
                  w[i,j] <- dat$w[dat_index]
                  w_ave[i,2] = w_ave[i,2] + 1
            }
      }
      u_ave[i,1] <- mean(u[i,])
      v_ave[i,1] <- mean(v[i,])
      w_ave[i,1] <- mean(w[i,])
      if (u_ave[i,2]==sampling_rate) {
            for (j in 1:sampling_rate) {
                  u_prime[i,j] <- u[i,j] - u_ave[i,1]
                  v_prime[i,j] <- v[i,j] - v_ave[i,1]
                  w_prime[i,j] <- w[i,j] - w_ave[i,1]
            }
      }
      if (u_ave[i,2]==sampling_rate) {
            uu[i,1] <- mean((u_prime[i,]^2))
            vv[i,1] <- mean((v_prime[i,]^2))
            ww[i,1] <- mean((w_prime[i,]^2))
            uv[i,1] <- mean((u_prime[i,]*v_prime[i,]))
            uw[i,1] <- mean((u_prime[i,]*w_prime[i,]))
            vw[i,1] <- mean((v_prime[i,]*w_prime[i,]))
      }
}

par(mfrow = c(3,1), mar = c(6,6,3,3))
# REM: x-axis is datetime range
plot(datetime[,2], u_ave[,1], ylim = c(-0.5, 0.5), xlim = c(77250, 77500), type = "l",ylab = "u (m/s)", xlab = "")
plot(datetime[,2], v_ave[,1], ylim = c(-0.5, 1), xlim = c(77250, 77500), type = "l",ylab = "v (m/s)", xlab = "")
plot(datetime[,2], w_ave[,1], ylim = c(-1, 1), xlim = c(77250, 77500), type = "l",ylab = "w (m/s)", xlab = "Time (s)")

# to zoom in on an area of interest, find indices:
start <- which(datetime[,2]==71450)
stop <- which(datetime[,2]==72000)
dt_zoom <- datetime[start:stop,2]
u_ave_zoom <- u_ave[start:stop,1]

par(mfrow = c(3,1), mar = c(4,4,2,2))
plot(datetime[,2], uu, ylim = c(0, 0.2), xlim = c(77250, 77500), type = "l", ylab = "uu", xlab = "")
plot(datetime[,2], vv, ylim = c(0, 1), xlim = c(77250, 77500), type = "l", ylab = "vv", xlab = "")
plot(datetime[,2], ww, ylim = c(0, 1), xlim = c(77250, 77500), type = "l", ylab = "ww", xlab = "Time (s)")

par(mfrow = c(3,1), mar = c(4,4,2,2))
plot(datetime[,2], uv, ylim = c(-0.3, 0.3), xlim = c(77250, 77500), type = "l", ylab = "uv", xlab = "")
plot(datetime[,2], uw, ylim = c(-0.3, 0.3), xlim = c(77250, 77500), type = "l", ylab = "uw", xlab = "")
plot(datetime[,2], vw, ylim = c(-0.3, 0.3), xlim = c(77250, 77500), type = "l", ylab = "vw", xlab = "Time (s)")

# Spectra
start <- match(77300,datetime[,2])
stop <- match(77340,datetime[,2])
par(mfrow = c(3,1), mar = c(4,4,2,2))
hist(uv[start:stop,1], breaks = c(-10000,-1.5,-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1.5,10000), xlim = c(-1.5,1.5), ylab = "u'v'", xlab = "", main = "")
hist(uw[start:stop,1], breaks = c(-10000,-1.5,-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1.5,10000), xlim = c(-1.5,1.5), ylab = "u'w'", xlab = "", main = "")
hist(vw[start:stop,1], breaks = c(-10000,-1.5,-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1.5,10000), xlim = c(-1.5,1.5), ylab = "v'w'", xlab = "", main = "")


