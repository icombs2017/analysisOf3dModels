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
ggsave("prism_error_means.pdf", plot= plot, width=12, height=4, units="in", dpi=300)