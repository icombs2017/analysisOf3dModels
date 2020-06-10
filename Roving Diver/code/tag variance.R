library(reshape2)
library(ggplot2)

#run t test between old and new cattle tag

tag <- read.csv("../data/cowtag areas.csv", head = T)
tag
tag$X8.24.18 <- tag$X8.24.18*10000
tag$X9.11.18 <- tag$X9.11.18*10000
tag$X11.8.19 <- tag$X11.8.19*10000

tag.melt=melt(tag,id.vars=c(1))
head(tag.melt)
colnames(tag.melt)[2]="date"
colnames(tag.melt)[3]="area"
head(tag.melt)

mean(tag.melt$area)
sd <- sd(tag.melt$area)
error <- sd/sqrt(84)
error

j <- read.csv("../data/imageJ.csv", head = T)

mean(j$area)
sd.j <- sd(j$area)
error.j <- sd.j/sqrt(9)
error.j
t.test(tag.melt$area, j$area)
