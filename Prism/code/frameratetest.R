# install.packages("ggpubr")
# install.packages("plyr")
# install.packages("Rmisc")
#library(plyr)
# Developmental version
#devtools::install_github("karthik/wesanderson")
# CRAN version
#install.packages("wesanderson")
library(ggplot2)
library(stringr)
library(gridExtra)
library(ggpubr)
library(Rmisc)
library(FSA)
library(rcompanion)
library(RColorBrewer)
library(dplyr)
library(vegan)
library(MANOVA.RM)
library(RVAideMemoire)
library(pairwiseAdonis)

frame <- read.csv("../data/frameratetest.new.csv", header = TRUE)

frame

# frame.melt <- melt(frame, id.vars=c(1:3), na.rm = T)
# head(frame.melt)
colnames(frame)[5]="total"
colnames(frame)[6]="cowtag"
head(frame)



##################################################################
 
#Anova on frame rate

frame.kw <- kruskal.test(total ~ frame.rate, data = frame)
  compare_means(total ~ frame.rate, data=frame, method = "kruskal.test")
frame.kw

cowtag.kw <- kruskal.test(cowtag ~ frame.rate, data = frame)
cowtag.kw

