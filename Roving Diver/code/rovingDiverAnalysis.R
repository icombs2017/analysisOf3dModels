
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, stringr, gridExtra, ggpubr, Rmisc, FSA, rcompanion, RColorBrewer, dplyr, vegan, nparcomp, RVAideMemoire, MANOVA.RM, pairwiseAdonis, PMCMR, PMCMRplus, patchwork)
pacman::p_load_gh("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
if(!require("devtools")) install.packages("devtools")
require(devtools)
install_github('davidcarslaw/openair')


# Making data the way i need it

rd <-read.csv("../data/RovingDiverSurveysFixed.csv", head = T, na.strings = "TMTC")
head(rd)
class(rd$TotalHealthy)

#Subset just for manuscript


rd$MonthYear = format(as.Date(as.character(rd$Date), format = "%m/%d/%y"), "%m-%y")
rd$Date <- as.Date(as.character(rd$Date), "%m/%d/%y")


rd.man <- subset(rd, Date < "2019-07-03") %>%  
  subset(SiteCode %in% c("SEFL01", "SEFL02", "SEFL03", "SEFL04", "SEFL05", "SEFL06", "SEFL08", "SEFL11", "SEFL12", "SLR North", "SLR South", "SLR Central", "SLR Ledge", "BC1", "T328", "FTL4")) 
  #subset(TotalHealthy != 'TMTC') 
head(rd.man)
class(rd.man$TotalHealthy)





#write.csv(rd.man, "rd.man.csv")

# summary(rd.man.1.1$TotalHealthy)
# head(rd.man.1.1)
# rd.man.1$SiteCode <- droplevels(rd.man.1$SiteCode)
# summary(rd.man.1$SiteCode)
# rd.man.1 <- rd.man.1.1
# head(rd.man.1)
# rd.man.1$prevalence <- (rd.man.1$TotalDiseased) / (rd.man.1$TotalDiseased + rd.man.1$TotalHealthy)
# head(rd.man.1)
# rd.man.1$Date2 <- NULL
# rd.man.1$Notes <- NULL
# head(rd.man.1)
# write.csv(rd.man.1, "rd.man.1.csv")
# rd.man <- read.csv("rd.man.1.csv", head = TRUE)
# head(rd.man)

# colnames(rd.man)[23]="Diseased"
# colnames(rd.man)[24]="Healthy"
# colnames(rd.man)[6]="County"
# rd.man$Diseased <- as.numeric(rd.man$Diseased)
# rd.man$Healthy <- as.numeric(rd.man$Healthy)
# rd.man$Total <- rd.man$Diseased + rd.man$Healthy
# head(rd.man)
# write.csv(rd.man, "rd.man.csv")
# 
# rd.man <- read.csv("rd.man.csv", head = TRUE)
# head(rd.man)
#rd.man <-read.csv("rd.man.csv", header = T, na.strings = 'TMTC')
rd.man$Total <- rd.man$UK + rd.man$DS + rd.man$BB + rd.man$RB + rd.man$YB + rd.man$TLD + rd.man$WP + rd.man$WS + rd.man$P + rd.man$PB + rd.man$BL + rd.man$TotalHealthy

head(rd.man)

countyData = rd.man %>% dplyr::select(c(County, MonthYear, SpeciesCode, TLD, Total, SiteCode)) %>% group_by(SiteCode,MonthYear,County) %>% summarise_if(is.integer, sum)


countyData$Prevalence = countyData$TLD/countyData$Total
head(countyData)


set.seed(999)
adonis(formula = Prevalence ~ County*MonthYear, data = countyData, method = "euclidian")

prev.pair <- pairwise.adonis(countyData[c(6)],factors=countyData$County,
                             sim.method='euclidian',p.adjust.m='bonferroni', perm = 999)

prev.pair



martin <- subset(countyData, County == "Martin")
palm <- subset(countyData, County == "Palm Beach")
broward <- subset(countyData, County == "Broward")

mean(martin$Total)
sd.martin <- sd(martin$Total)
len.martin <-length(martin$Total)
error.martin <- sd.martin/sqrt(len.martin)
error.martin

mean(palm$Total)
sd.palm <- sd(palm$Total)
len.palm <-length(palm$Total)
error.palm <- sd.palm/sqrt(len.palm)
error.palm

mean(broward$Total)
sd.broward <- sd(broward$Total)
len.broward <-length(broward$Total)
error.broward <- sd.broward/sqrt(len.broward)
error.broward

