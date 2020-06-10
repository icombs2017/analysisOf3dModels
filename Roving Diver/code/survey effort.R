
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, stringr, gridExtra, ggpubr, Rmisc, FSA, rcompanion, RColorBrewer, dplyr, vegan, nparcomp, RVAideMemoire, MANOVA.RM, pairwiseAdonis, PMCMR, PMCMRplus, patchwork, reshape2)
pacman::p_load_gh("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

rd.old <- read.csv("../data/survey.effort.csv", head = T)
head(rd.old)



#melt data so all numbers are in one column
rd.melt <- melt(rd.old,id.vars=c(1:3), na.rm = TRUE)
head(rd.melt)

colnames(rd.melt)[4]="date"
colnames(rd.melt)[5]="prevalence"
head(rd.melt)


#Checking Normality
shapiro.test(rd.melt$prevalence)

# Shapiro-Wilk normality test
# 
# data:  rd.melt$prevalence
# W = 0.63313, p-value = 8.857e-11


#Non-normal so going with univariate PERMANOVA
#Design is not conducive to Friedman's test
set.seed(999)
adonis(formula = prevalence ~ location*date, data = rd.melt, method = "euclidian")

prev.pair <- pairwise.adonis(rd.melt[c(5)],factors=rd.melt$county,
               sim.method='euclidian',p.adjust.m='bonferroni', perm = 999)

pair.prev <- rbind(prev.pair)

write.csv(pair.prev, file="roving.diver.pairwise.csv")




aov.melt <- aov(prevalence ~ county*date, data = rd.melt)
summary(aov.melt)
# Going to confirm with a KW

prev.kw.site <- compare_means(prevalence ~ date, data=rd.melt, method = "kruskal.test")
prev.kw.site
prev.dunn.site <- dunnTest(prevalence ~ date, data=rd.melt, method = "bh")$res
prev.dunn.site




prev.kw.date <- compare_means(prevalence ~ date, data=rd.melt, method = "kruskal.test")
prev.kw.date
prev.dunn.date <- dunnTest(prevalence ~ date, data=rd.melt, method = "bh")$res
prev.dunn.date




prev.kw.county <- compare_means(prevalence ~ county, data=rd.melt, method = "kruskal.test")
prev.kw.county
prev.dunn.county <- dunnTest(prevalence ~ county, data=rd.melt, method = "bh")$res
prev.dunn.county







########################################################################
########################################################################
#########################################################################
#########################################################################
##########################################################################
###########################################################################
###########################################################################


# Subset for county
martin <- subset(rd.melt, county == "Martin")
#martin <- subset(martin.1, date != "Nov.17")
droplevels(martin$county)
martin
palm.beach <- subset(rd.melt, county == "Palm Beach")
droplevels(palm.beach$county)
palm.beach
broward <- subset(rd.melt, county == "Broward")
droplevels(broward$county)
broward



ggplot(rd.melt, aes(x=date, y=prevalence, color = county, group = county, label = county))+
  geom_line()+
  geom_point()
  
# geom_label()+
  # geom_text()+
  #geom_errorbar()

martin.prev.plot <- 
  ggplot(martin, aes(x=date, y=prevalence, color = site, group = site, label=site))+
  geom_line(size = 1)+
  geom_point(size = 2)+
  #stat_summary(geom = 'errorbar', fun.data = mean_se, position = "dodge")+
  # geom_label()+
  # geom_text()+
  labs(y = "Prevalence")+
  scale_x_discrete(limits=c("Nov.17","Dec.17","Jan.18", "Apr.18", "May.18", "Jun.18", "Aug.18", "Nov.18", "Dec.18", "Feb.19", "Mar.19", "Apr.19", "May.19", "Jun.19"), labels =c("Nov 17","Dec 17","Jan 18", "Apr 18", "May 18", "Jun 18", "Aug 18", "Nov 18", "Dec 18", "Feb 19", "Mar 19", "Apr 19", "May 19", "Jun 19"))+
  ylim(0,0.5)

# takes the legend and saves it as a separate object (grob)
get_legend<-function(martin.prev.plot){
  tmp <- ggplot_gtable(ggplot_build(martin.prev.plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  martin.legend <- tmp$grobs[[leg]]
  return(martin.legend)
}
martin.legend=get_legend(martin.prev.plot)


#theme(legend.position = "none")

martin.prev.plot <- martin.prev.plot + theme(
  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),                 panel.background = element_rect(fill = '#F5F5F5'),
  plot.title = element_text(hjust = 0.5), 
  axis.line = element_line(colour = "black"), 
  axis.ticks = element_line(color="black"), 
  text = element_text(size=15, color="black"), 
  # axis.text.x=element_text(angle = 45, hjust = 1,size=12, color="black"), 
  # axis.text.y=element_text(size=12, color="black")
  axis.title.x=element_blank(),
  axis.text.x=element_blank(#angle = 45, hjust = 1,size=12, color="black"),
    ),
  axis.text.y=element_text(size=12, color="black")
  )
martin.prev.plot <- martin.prev.plot + rremove("legend")









palm.prev.plot <- 
  ggplot(palm.beach, aes(x=date, y=prevalence, color = site, group = site, label = site))+
  geom_line()+
  geom_point()+
  # geom_label()+
  # geom_text()+
  scale_x_discrete(limits=c("Nov.17","Dec.17","Jan.18", "Apr.18", "May.18", "Jun.18", "Aug.18", "Nov.18", "Dec.18", "Feb.19", "Mar.19", "Apr.19", "May.19", "Jun.19"), labels =c("Nov 17","Dec 17","Jan 18", "Apr 18", "May 18", "Jun 18", "Aug 18", "Nov 18", "Dec 18", "Feb 19", "Mar 19", "Apr 19", "May 19", "Jun 19"))+
  ylim(0,0.5)+
  labs(y = "Prevalence", x = "Month")

# takes the legend and saves it as a separate object (grob)
get_legend<-function(palm.prev.plot){
  tmp <- ggplot_gtable(ggplot_build(palm.prev.plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  palm.legend <- tmp$grobs[[leg]]
  return(palm.legend)
}
palm.legend=get_legend(palm.prev.plot)


#theme(legend.position = "none")

palm.prev.plot <- palm.prev.plot+theme(
  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),                 panel.background = element_rect(fill = '#F5F5F5'),
  plot.title = element_text(hjust = 0.5), 
  axis.line = element_line(colour = "black"), 
  axis.ticks = element_line(color="black"), 
  text = element_text(size=15, color="black"), 
  # axis.text.x=element_text(angle = 45, hjust = 1,size=12, color="black"), 
  # axis.text.y=element_text(size=12, color="black")
  axis.title.x=element_blank(),
  axis.text.x=element_blank(#angle = 45, hjust = 1,size=12, color="black"),
  ),
  axis.text.y=element_text(size=12, color="black")
)
palm.prev.plot <- palm.prev.plot + rremove("legend")





broward.prev.plot <- 
  ggplot(broward, aes(x=date, y=prevalence, color = site, group = site, label = site))+
  geom_line()+
  geom_point()+
  # geom_label()+
  # geom_text()+
  scale_x_discrete(name = "Month",limits=c("Nov.17","Dec.17","Jan.18", "Apr.18", "May.18", "Jun.18", "Aug.18", "Nov.18", "Dec.18", "Feb.19", "Mar.19", "Apr.19", "May.19", "Jun.19"), labels =c("Nov 17","Dec 17","Jan 18", "Apr 18", "May 18", "Jun 18", "Aug 18", "Nov 18", "Dec 18", "Feb 19", "Mar 19", "Apr 19", "May 19", "Jun 19"))+
  ylim(0,0.5)+
  annotate("label", x = 2.45, y = 0.45, label = paste ("Broward County"), size = 5) +
  labs(y = "Prevalence", x = "Month")

 
get_legend<-function(broward.prev.plot){
  tmp <- ggplot_gtable(ggplot_build(broward.prev.plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  broward.legend <- tmp$grobs[[leg]]
  return(broward.legend)
}
broward.legend=get_legend(broward.prev.plot)


 #theme(legend.position = "none")

broward.prev.plot <- broward.prev.plot+theme(
  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),                 panel.background = element_rect(fill = '#F5F5F5'),
  plot.title = element_text(hjust = 0.5), 
  axis.line = element_line(colour = "black"), 
  axis.ticks = element_line(color="black"), 
  text = element_text(size=15, color="black"), 
  axis.text.x=element_text(angle = 45, hjust = 1,size=12, color="black"), 
  axis.text.y=element_text(size=12, color="black")
                    )

broward.prev.plot <- broward.prev.plot + rremove("legend")

prev.multi=grid.arrange(martin.prev.plot, martin.legend,palm.prev.plot, palm.legend, broward.prev.plot,broward.legend, ncol=2, nrow=3, widths=c(3,0.75), heights=c(3,3,4.1))

ggsave("site.prevalence.eps", plot= prev.multi, width=12, height=8, units="in", dpi=600)

















ggboxplot(
  rd.melt,
  x = "county",
  y = "prevalence",
  #color = "grey30",
  #palette = c("#7BA46C", "#EACF9E", "#008D91"),
  #fill = "site",
  add = "jitter",
  add.params = list(size = 1, jitter = 0.5),
  width = 0.7,
  size = 0.5
) + labs(x = "County",
         y = "Average Prevalence",
         fill = 'County') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none") + 
  scale_x_discrete(limits=c("Broward", "Palm Beach", "Martin"))+
  ylim(0,0.5)+theme(plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black"), axis.ticks = element_line(color="black"), text = element_text(size=15, color="black"), axis.text.x=element_text(size=12, color="black"), axis.text.y=element_text(size=12, color="black")) +
  stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=2,hjust=0.5, vjust=4)
  #geom_text(data=lesion.count.list, aes(x = Group, y=0, vjust=-22.5, label=Letter)) +
  #scale_y_continuous(labels = scale)



ggplot(rd.melt, aes(x=county, y=prevalence, color = county, group = county, label = county))+
  geom_boxplot()+
  geom_point()+
  scale_x_discrete(limits=c("Martin", "Palm Beach", "Broward"))+
  ylim(0,0.5)+
  theme(plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black"), axis.ticks = element_line(color="black"), text = element_text(size=15, color="black"), axis.text.x=element_text(size=12, color="black"), axis.text.y=element_text(size=12, color="black")) +
  labs(x = "County",
       y = "Average Disease Prevalence",
       fill = 'County') 
# geom_label()+
# geom_text()+
#geom_errorbar()


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################


####### Now running pevalence numbers by county and by month

prev <- read.csv("prev.county.csv", header = TRUE)
head(prev)


error.prev <- ggplot(prev, aes(date, prevalence, group = county, fill = county, color = county))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin = prevalence - stderror, 
                    ymax = prevalence + stderror),
                width = 0.5)+
  #stat_summary(geom = 'errorbar', fun.data = mean_se, position = "dodge")+
  scale_x_discrete(name = "Month",limits=c("Nov.17","Dec.17","Jan.18", "Apr.18", "May.18", "Jun.18", "Aug.18", "Nov.18", "Dec.18", "Feb.19", "Mar.19", "Apr.19", "May.19", "Jun.19"), labels =c("Nov 17","Dec 17","Jan 18", "Apr 18", "May 18", "Jun 18", "Aug 18", "Nov 18", "Dec 18", "Feb 19", "Mar 19", "Apr 19", "May 19", "Jun 19"))+
  ylim(0,0.5)+
  labs(y = "Prevalence", x = "Month")
error.prev

error.prev <- error.prev+ theme(
  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),                 panel.background = element_rect(fill = '#F5F5F5'), 
  plot.title = element_text(hjust = 0.5), 
  axis.line = element_line(colour = "black"), 
  axis.ticks = element_line(color="black"), 
  text = element_text(size=20, color="black"), 
  axis.text.x=element_text(angle = 45, hjust = 1, size=14, color="black"), 
  axis.text.y=element_text(size=12, color="black")
)
ggsave("prev.county.eps", plot= error.prev, width=12, height=8, units="in", dpi=600)


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################



####### Now looking at global trends across NRFT


g.prev <- read.csv("global.prev.csv", header = TRUE)
head(g.prev)

shapiro.test(g.prev$prevalence)



g.error.prev <- ggplot(g.prev, aes(date, prevalence, group = location, fill = location, color = location))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin = prevalence - stderror, 
                    ymax = prevalence + stderror), 
                width = .5)+
  theme(plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black"), axis.ticks = element_line(color="black"), text = element_text(size=15, color="black"), axis.text.x=element_text(size=12, color="black"), axis.text.y=element_text(size=12, color="black")) +
  #stat_summary(geom = 'errorbar', fun.data = mean_se, position = "dodge")+
  scale_x_discrete(name = "Month",limits=c("Nov.17","Dec.17","Jan.18", "Apr.18", "May.18", "Jun.18", "Aug.18", "Nov.18", "Dec.18", "Feb.19", "Mar.19", "Apr.19", "May.19", "Jun.19"), labels =c("Nov 17","Dec 17","Jan 18", "Apr 18", "May 18", "Jun 18", "Aug 18", "Nov 18", "Dec 18", "Feb 19", "Mar 19", "Apr 19", "May 19", "Jun 19"))+
  ylim(0,0.5)+
  labs(y = "Prevalence", x = "Month")
g.error.prev

g.error.prev <- g.error.prev + theme(
                panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),                 panel.background = element_rect(fill = '#F5F5F5'), 
                plot.title = element_text(hjust = 0.5), 
                axis.line = element_line(colour = "black"), 
                axis.ticks = element_line(color="black"), 
                text = element_text(size=20, color="black"), 
                axis.text.x=element_text(angle = 45, hjust = 1, size=14, color="black"), 
                axis.text.y=element_text(size=12, color="black")
                )




ggsave("prev.global.eps", plot= g.error.prev, width=12, height=8, units="in", dpi=600)

#gainsboro = #DCDCDC	
#light gray = #D3D3D3	
#white smoke = #F5F5F5




martin.prev.plot <- 
  ggplot(martin, aes(x=date, y=prevalence, color = site, group = site, label=site))+
  geom_line(size = 1)+
  geom_point(size = 2)+
  #stat_summary(geom = 'errorbar', fun.data = mean_se, position = "dodge")+
  # geom_label()+
  # geom_text()+
  labs(y = "Prevalence")+
  scale_x_discrete(limits=c("Mar.19", "Apr.19", "May.19", "Jun.19"), labels =c("Mar 19", "Apr 19", "May 19", "Jun 19"))+
  ylim(0,0.5)

ggsave("slr_permit.eps", plot = martin.prev.plot, width = 12, height = 8, units = 'in', dpi = 600)


theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
      panel.background = element_rect(fill = '#F5F5F5'),
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(color="black"),
      text = element_text(size=15, color="black"),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=12, color="black")
)




martin.prev.plot <- martin.prev.plot1 + 
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"), 
        panel.background = element_rect(fill = '#F5F5F5'),
        plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = "black"), 
        axis.ticks = element_line(color="black"), 
        text = element_text(size=15, color="black"), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=12, color="black")
  )







