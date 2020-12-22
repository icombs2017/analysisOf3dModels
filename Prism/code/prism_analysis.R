####packages####

library(vegan)
library(parallel)
library(pairwiseAdonis)
library(ggpubr)
library(FSA)
library(dplyr)
library(rcompanion)
library(ggplot2)
library(gridExtra)


####data import####

# if returning to previously-run stats
load("prism_stats.RData")

prism <- read.csv("prism_area_error.csv", head=T)
head(prism)
str(prism)

# converting sq m to sq cm
prism$square <- prism$square*10000
prism$rectangle <- prism$rectangle*10000
prism$circle <- prism$circle*10000
head(prism)

# calculating absolute difference between measured shape areas and the actual template
prism$square.diff <- abs(prism$square-40.3225)
prism$rectangle.diff <- abs(prism$rectangle-12.9032)
prism$circle.diff <- abs(prism$circle-45.60367312)
head(prism)

# calculating error between measured shape areas and the actual template
prism$square.error <- abs((prism$square-40.3225)/40.3225)
prism$rectangle.error <- abs((prism$rectangle-12.9032)/12.9032)
prism$circle.error <- abs((prism$circle-45.60367312)/45.60367312)
head(prism)


####shape error normality####

square.norm <- shapiro.test(prism$square.diff)
square.norm

rectangle.norm <- shapiro.test(prism$rectangle.diff)
rectangle.norm

circle.norm <- shapiro.test(prism$circle.diff)
circle.norm

shapiro <- rbind(square.norm,rectangle.norm,circle.norm)
write.csv(shapiro, file="Shapiro_error_outputs.csv")

# histograms
square.hist <- histogram(prism$square.diff)
rectangle.hist <- histogram(prism$rectangle.diff)
circle.hist <- histogram(prism$circle.diff)
# very right-skewed, applying a sqrt transformation

# very non-normal, appyling sqrt transformation
prism$square.diff.sqrt <- sqrt(prism$square.diff)
prism$rectangle.diff.sqrt <- sqrt(prism$rectangle.diff)
prism$circle.diff.sqrt <- sqrt(prism$circle.diff)
head(prism)
str(prism)

# histograms
square.sqrt.hist <- histogram(prism$square.diff.sqrt)
rectangle.sqrt.hist <- histogram(prism$rectangle.diff.sqrt)
circle.sqrt.hist <- histogram(prism$circle.diff.sqrt)
# less right-skew

hist = grid.arrange(square.hist, rectangle.hist, circle.hist, square.sqrt.hist, rectangle.sqrt.hist, circle.sqrt.hist, ncol=3, nrow=2, widths=c(3,3,3), heights=c(3,3))
ggsave("prism_error_histogram.pdf", plot= hist, width=9, height=6, units="in", dpi=300)

# Shapiro test of data normality
square.sqrt.norm <- shapiro.test(prism$square.diff.sqrt)
square.sqrt.norm

rectangle.sqrt.norm <- shapiro.test(prism$rectangle.diff.sqrt)
rectangle.sqrt.norm

circle.sqrt.norm <- shapiro.test(prism$circle.diff.sqrt)
circle.sqrt.norm

shapiro.sqrt <- rbind(square.sqrt.norm,rectangle.sqrt.norm,circle.sqrt.norm)
write.csv(shapiro.sqrt, file="Shapiro_error_sqrt_outputs.csv")
# still not normal, proceeding with nonparametric tests


####shape error PERMANOVA####

# Setting seed allows randomized processes to be repeated later
set.seed(215)

# creates a dissimilarity matrix of square-root transformed shape areas, based on Euclidean distance
error.dist <- vegdist(prism[c(14:16)], method="euclidean", na.rm = TRUE)
# replaces missing data (NAs) with 0
error.dist[is.na(error.dist)] <- 0

# testing for homogeneity of variance among sites  (pool vs BC1 vs T328 vs FTL4)
error.disp <- betadisper(error.dist, group=prism$site)
error.site.disp <- permutest(error.disp, bias.adjust = TRUE, perm = 9999)
error.site.disp
# copy output to a csv file

# post hoc tests for site
error.site.disp.HSD <- TukeyHSD(error.disp)
error.site.disp.HSD
# indicates unequal variance between pool vs FTL4 and T328 vs FTL4
# copy output to a csv file

# plots showing  variance among sites
plot(error.disp, ellipse = TRUE, hull = FALSE)
boxplot(error.disp)
# save both plots

# running the PERMANOVA
error.perm <- adonis(error.dist ~ site, data=prism, permutations = 1e6, parallel = getOption("mc.cores"))
error.perm
# No significant differences in shape errors among sites

# exporting PERMANOVA results
write.csv(error.perm$aov.tab, file = "PERMANOVA_error_outputs.csv")


####exporting/importing results####

# saving dataframes and test results
save(prism, error.dist, error.perm, file = "prism_stats.RData")
# to load back in later if needed
load("prism_stats.RData")


####shape error plots####

prism$site=factor(prism$site, levels=unique(prism$site)) 

# creates a function that forces y axis labels to have 1 decimal place
scale <- function(x) sprintf("%.0f", x)

# boxplots comparing shape errors across sites
square.diff.box <-
  ggboxplot(
    prism,
    x = "site",
    y = "square.diff",
    color = "grey30",
    palette = c("grey", "#7BA46C", "#EACF9E", "#008D91"),
    fill = "site",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Site",
           y = "Shape Error (cm2)",
           title = "Square",
           fill = 'Site') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none") + 
  scale_y_continuous(labels = scale)
square.diff.box

rectangle.diff.box <-
  ggboxplot(
    prism,
    x = "site",
    y = "rectangle.diff",
    color = "grey30",
    palette = c("grey", "#7BA46C", "#EACF9E", "#008D91"),
    fill = "site",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Site",
           y = element_blank(),
           title = "Rectangle",
           fill = 'Site') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none") + 
   scale_y_continuous(labels = scale)
rectangle.diff.box

circle.diff.box <-
  ggboxplot(
    prism,
    x = "site",
    y = "circle.diff",
    color = "grey30",
    palette = c("grey", "#7BA46C", "#EACF9E", "#008D91"),
    fill = "site",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Site",
           y = element_blank(),
           title = "Circle",
           fill = 'Site') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none") + 
   scale_y_continuous(labels = scale)
circle.diff.box

plot=grid.arrange(square.diff.box, rectangle.diff.box, circle.diff.box, ncol=3, nrow=1, widths=c(3.1,3,3), heights=c(3))

#saves plot as PDF
<<<<<<< HEAD
ggsave("prism_error_means.pdf", plot= plot, width=12, height=4, units="in", dpi=300)
=======
ggsave("prism_error_means.pdf", plot= plot2, width=12, height=4, units="in", dpi=300)

#-------------------------------
####shape error correlations####

# testing correlation between size of model (total area) and error associated with shape measurements
# i.e. if the model area is larger (higher filming altitude), are models worse quality?

# tests correlation between model area and the difference between in-model shape measurements and the template measurements (error)
cor.test(prism.sub$whole,prism.sub$square.diff, method="spearman", use="complete.obs")
# Spearman's rank correlation rho
# 
# data:  prism.sub$whole and prism.sub$square.diff
# S = 3958750, p-value = 0.1084
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.09333898 

cor.test(prism.sub$whole,prism.sub$rectangle.diff, method="spearman", use="complete.obs")
# Spearman's rank correlation rho
# 
# data:  prism.sub$whole and prism.sub$rectangle.diff
# S = 3100600, p-value = 1.231e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2527713 

cor.test(prism.sub$whole,prism.sub$circle.diff, method="spearman", use="complete.obs")
# Spearman's rank correlation rho
# 
# data:  prism.sub$whole and prism.sub$circle.diff
# S = 48656, p-value = 0.01591
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2794314 

# significant correlation for rectangle and circle, so let's see if this is driven by particular sites

#-----------------------------------
####shape error correlation plots####

# correlation plots comparing shape error to overall model area

prism.sub$location=factor(prism.sub$location, levels=unique(prism.sub$location)) 

# makes a scatterplot of correlation
corr.square <-
  ggscatter(
    prism.sub,
    x = "whole",
    y = "square.diff",
    color = "location",
    add = "reg.line",
    conf.int = TRUE,
    palette = c("grey", "#3070cf", "#79d7fb"),
    label.x.npc = "middle",
    label.y.npc = "top",
    xlab = "Model Area (m2)",
    ylab = "Shape Error (cm2)"
    #main = "Square"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~','~"), color=location))
corr.square

# takes the legend and saves it as a separate object (grob)
get_legend<-function(corr.square){
  tmp <- ggplot_gtable(ggplot_build(corr.square))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend=get_legend(corr.square)

# removes legend
corr.square <- corr.square + theme_bw()
corr.square <- corr.square + rremove("legend")

corr.rectangle <-
  ggscatter(
    prism.sub,
    x = "whole",
    y = "rectangle.diff",
    color = "location",
    add = "reg.line",
    conf.int = TRUE,
    palette = c("grey", "#3070cf", "#79d7fb")
    #main = "Rectangle"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~','~"), color=location)) + 
  labs(x = "Model Area (m2)",
       y = element_blank(),
       #title = "Circle",
       fill = 'Location') +
  theme_bw() +
  rremove("legend")
corr.rectangle

corr.circle <-
  ggscatter(
    prism.sub,
    x = "whole",
    y = "circle.diff",
    color = "location",
    add = "reg.line",
    conf.int = TRUE,
    palette = c("grey", "#3070cf", "#79d7fb"),
    xlab = "Model Area (m2)",
    ylab = FALSE
    #main = "Circle"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~','~"), color=location)) + 
  labs(x = "Model Area (m2)",
       y = element_blank(),
       #title = "Circle",
       fill = 'Location') +
  theme_bw() + 
  rremove("legend")
corr.circle

plot3=grid.arrange(corr.square, corr.rectangle, corr.circle, legend, ncol=3, nrow=2, layout_matrix=rbind(c(1,2,3),c(4,4,4)), widths=c(3.1,3,3), heights=c(3,0.25))

#saves plot as PDF (transparent CI won't export as EPS)
ggsave("corr_area_error.pdf", plot= plot3, width=12, height=4, units="in", dpi=300)

# so, there may be a sig increase in shape measurement error as model area increases, but depending on location (stlucie)

#-----------------------------------
####all figure multiplot####

# multiplot of all generated figures

multiplot=grid.arrange(square.box, rectangle.box, circle.box, square.diff.box, rectangle.diff.box, circle.diff.box, corr.square, corr.rectangle, corr.circle, legend, ncol=3, nrow=4, layout_matrix=rbind(c(1,2,3),c(4,5,6),c(7,8,9),c(10,10,10)), widths=c(3.1,3,3), heights=c(3,3,3,0.25))

ggsave("prism_area_error_panel.pdf", plot= multiplot, width=12, height=12, units="in", dpi=300)

#-----------------------------------

########## Testing for differences in model accuracy based on different frame rates

## A subset of model were extracted across their time points at 3 different framerates: 3, 4 and 5 fps


frame <- read.csv("../data/frameratetest.csv", header = TRUE)

frame

# Kruskal-Wallis on frame rate

frame.kw <- kruskal.test(total ~ frame.rate, data = frame)
frame.kw

# no differences in model area based on framerate!


>>>>>>> 2457e76834eba310a67146b6ef5ea37c0e14449c
