rd <- read.xls("../data/RovingDiverSurveysManuscript.xlsx", head = T, na.strings ="TMTC")
class(rd$Healthy)
head(rd)


rd$MonthYear = format(as.Date(as.Date(rd$Date), format = "%y/%m/%d"), "%m-%y")
rd$Date = as.Date(rd$Date)
class(rd$Date)
head(rd)


rd$County = if_else(rd$SiteCode %in% c("SEFL01", "SEFL02", "SEFL03", "SLR South", "SLR Central", "SLR Ledge", "SLR North"), "Martin",
                    if_else(rd$SiteCode %in% c("T328", "BC1", "FTL4"), "Broward", "Palm Beach"))

rd$Location= if_else(rd$SiteCode %in% c("SEFL01", "SEFL02", "SEFL03", "SLR South", "SLR Central", "SLR Ledge", "SLR North"), "SLR",
                     if_else(rd$SiteCode %in% c("T328", "BC1", "FTL4"), "PMP",
                             if_else(rd$SiteCode %in% c("SEFL04", "SEFL05", "SEFL06"), "JUP",
                                     if_else(rd$SiteCode %in% c("SEFL08", "SEFL09", "SEFL10", "SEFL11", "SEFL12"), "WPB", ""))))
head(rd)



rd$TotalObservations <- rd$UK + rd$DS + rd$BB + rd$RB + rd$YB + rd$SCTLD + rd$WP + rd$WS + rd$P + rd$PB + rd$BL + rd$Healthy
head(rd)

rd.site <- subset(rd, SiteCode %in% c("SEFL01", "SEFL02", "SEFL04", "SEFL05", "SEFL06", "SEFL08", "SEFL11", "SEFL12", "SLR North", "SLR South", "SLR Central", "SLR Ledge", "BC1", "T328", "FTL4")) 

head(rd.site)
tail(rd.site)

siteData <- rd.site %>% dplyr::select(c(SiteCode, Location, County, MonthYear, SpeciesCode, SCTLD, TotalObservations)) %>%  
  group_by(SiteCode, MonthYear) %>% 
  summarise_if(is.integer, sum)
siteData

siteData$Prevalence = siteData$SCTLD/siteData$TotalObservations
head(siteData)

prev.dist <- vegdist(siteData$Prevalence, method="euclidean", na.rm= TRUE)
# replaces missing data (NAs) with 0
prev.dist[is.na(prev.dist)] <- 0

# testing for homogeneity of variance among locations
prev.disp <- betadisper(prev.dist, group=siteData$SiteCode)
permutest(prev.disp, bias.adjust = TRUE, perm = 9999)
# significant test value indicates heterogeneous variance

# post hoc tests
prev.disp.HSD <- TukeyHSD(prev.disp)
prev.disp.HSD
# indicates unequal variance in site SLR

# plots showing higher variance in Martin versus other locations
boxplot(prev.disp)
# These results suggest unequal variance among locations, with highest variance attributed to Martin
# Since Martin is also the smallest group, this can result in liberal PERMANOVA results (per Anderson and Walsh 2013)
# However, given the results of the PERMANOVA below, the interpretation remains that only models differ within locations

# running the PERMANOVA
set.seed(999)
prev.perm <- adonis(formula = Prevalence ~ SiteCode*MonthYear, data = siteData, method = "euclidian", permutations = 9999, na.rm = TRUE)
prev.perm









