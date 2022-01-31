#To analyze data for oxygen consumption equation

library(dplyr)
library(ggplot2)

setwd("/Users/alannabachtlin/Desktop/data/")
DO_BB629 <- data.frame(read.csv("DOBBlogprofile.csv"))

vel.dif <- array(NA, dim = c((nrow(DO_BB629)-1),(nrow(DO_BB629))))
do.dif <- vel.dif
height <- vel.dif
#u2-u1
for(j in 1:(nrow(DO_BB629)-1))  {
  for(i in (j+1):nrow(DO_BB629)) {
    vel.dif[j,i] <- DO_BB629$U.m.s.[i] - DO_BB629$U.m.s.[j]
    do.dif[j,i] <- DO_BB629$Avg.DO..mg.l.[i] - DO_BB629$Avg.DO..mg.l.[j]
    height[j,i] <- log(DO_BB629$z..m.[i] / DO_BB629$z..m.[j])
  }
}

plot(do.dif, height)
plot(vel.dif, height)

ggplot(DO_BB629) +
  geom_point(aes(x=Avg.DO..mg.l.,y=z..m.))
