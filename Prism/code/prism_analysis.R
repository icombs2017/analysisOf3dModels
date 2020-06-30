####packages####

library(ggplot2)
library(vegan)
library(gridExtra)
library(ggpubr)
library(FSA)
library(rcompanion)
library(vegan)
library(parallel)
library(pairwiseAdonis)

#------------------------------------
####prism model area####

model <- read.csv("prism_model_area.csv", head = T)
str(model)
head(model)

shapiro.test(model$area)
# Shapiro-Wilk normality test
# 
# data:  model$area
# W = 0.9114, p-value = 0.07856

# ANOVA of whole model areas

# testing for difference between model areas by location
model.anova <- aov(area~location, data = model)
summary(model.anova)

capture.output(summary(model.anova), file = "ANOVA_model_outputs.csv")

# calculates mean and SD by location
aggregate(area ~ location, model, mean)
# location     area
# 1     LBTS 8.934636
# 2     Pool 4.758737
# 3      SLR 5.500271
aggregate(area ~ location, model, sd)
# location      area
# 1     LBTS 1.2870938
# 2     Pool 0.8023665
# 3      SLR 1.7716244

# There is a significant difference between model areas, with reef sites having larger models, but larger standard deviation

#-----------------------------
####prism shape area####

prism <- read.csv("prism_area_error.csv", head=T)
head(prism)
str(prism)

# checking normality assumptions
shapiro.test(prism$square)
# Shapiro-Wilk normality test
# 
# data:  prism$square
# W = 0.91362, p-value = 4.538e-12
shapiro.test(prism$rectangle)
# Shapiro-Wilk normality test
# 
# data:  prism$rectangle
# W = 0.84311, p-value < 2.2e-16
shapiro.test(prism$circle)
# Shapiro-Wilk normality test
# 
# data:  prism$circle
# W = 0.90575, p-value = 3.736e-05

# histograms
hist(prism$square, breaks = 50)
hist(prism$rectangle, breaks = 50)
hist(prism$circle, breaks = 50)

# not normal, doing sqrt transformation due to right skew
prism$square.sqrt <- sqrt(prism$square)
prism$rectangle.sqrt <- sqrt(prism$rectangle)
prism$circle.sqrt <- sqrt(prism$circle)
head(prism)

shapiro.test(prism$square.sqrt)
# Shapiro-Wilk normality test
# 
# data:  prism$square.sqrt
# W = 0.92871, p-value = 9.45e-11
shapiro.test(prism$rectangle.sqrt)
# Shapiro-Wilk normality test
# 
# data:  prism$rectangle.sqrt
# W = 0.87949, p-value = 2.006e-14
shapiro.test(prism$circle.sqrt)
# Shapiro-Wilk normality test
# 
# data:  prism$circle.sqrt
# W = 0.91458, p-value = 8.984e-05

#------------------------------------
####shape area PERMANOVA####

# Setting seed allows randomized processes to be repeated later
set.seed(215)

# creates a dissimilarity matrix of square-root transformed shape areas, based on Euclidean distance
area.dist <- vegdist(prism[c(14:16)], method="euclidean", na.rm= TRUE)
# replaces missing data (NAs) with 0
area.dist[is.na(area.dist)] <- 0

# testing for homogeneity of variance among locations
area.disp <- betadisper(area.dist, group=prism$location)
permutest(area.disp, bias.adjust = TRUE, perm = 9999)
# significant test value indicates heterogeneous variance

# post hoc tests
area.disp.HSD <- TukeyHSD(area.disp)
area.disp.HSD
# indicates unequal variance in site SLR

# plots showing higher variance in SLR versus other locations
plot(area.disp, ellipse = TRUE, hull = FALSE)
boxplot(area.disp)
# These results suggest unequal variance among locations, with highest variance attributed to SLR
# Since SLR is also the smallest group, this can result in liberal PERMANOVA results (per Anderson and Walsh 2013)
# However, given the results of the PERMANOVA below, the interpretation remains that only models differ within locations

# running the PERMANOVA
# model is nested within location using '/', and is accounted as a random factor using 'strata'
# typically 9999 permutations are sufficient, but model terms were close to 0.05, so running with a million permutations
area.perm <- adonis(area.dist ~ location*location/model*face, data=prism, strata = prism$model, permutations = 1e6, parallel = getOption("mc.cores"))
area.perm

write.csv(area.perm$aov.tab, file = "PERMANOVA_area_outputs.csv")
# Significant differences between location and prism face, but not among models within location

# pairwise PERMANOVA based on significant factors (model)
model.pair <- pairwise.adonis2(area.dist ~ model, data = prism, strata = 'model', nperm = 1e6)
model.pair

model.pair.out <- bind_rows(model.pair$pool1_vs_pool3, model.pair$pool1_vs_pool4, model.pair$pool1_vs_pool5, model.pair$pool1_vs_pool6, model.pair$pool1_vs_pool7, model.pair$pool3_vs_pool4, model.pair$pool3_vs_pool5, model.pair$pool3_vs_pool6, model.pair$pool3_vs_pool7, model.pair$pool4_vs_pool5, model.pair$pool4_vs_pool6, model.pair$pool4_vs_pool7, model.pair$pool5_vs_pool6, model.pair$pool5_vs_pool7, model.pair$pool6_vs_pool7, model.pair$lbts1_vs_lbts2, model.pair$lbts1_vs_lbts3, model.pair$lbts1_vs_lbts5, model.pair$lbts1_vs_lbts6, model.pair$lbts1_vs_lbts7, model.pair$lbts1_vs_lbts8, model.pair$lbts2_vs_lbts3, model.pair$lbts2_vs_lbts5, model.pair$lbts2_vs_lbts6, model.pair$lbts2_vs_lbts7, model.pair$lbts2_vs_lbts8, model.pair$lbts3_vs_lbts5, model.pair$lbts3_vs_lbts6, model.pair$lbts3_vs_lbts7, model.pair$lbts3_vs_lbts8, model.pair$lbts5_vs_lbts6, model.pair$lbts5_vs_lbts7, model.pair$lbts5_vs_lbts8, model.pair$lbts6_vs_lbts7, model.pair$lbts6_vs_lbts8, model.pair$lbts7_vs_lbts8, model.pair$slr1_vs_slr2, model.pair$slr1_vs_slr6, model.pair$slr1_vs_slr7, model.pair$slr1_vs_slr8, model.pair$slr1_vs_slr9, model.pair$slr2_vs_slr6, model.pair$slr2_vs_slr7, model.pair$slr2_vs_slr8, model.pair$slr2_vs_slr9, model.pair$slr6_vs_slr7, model.pair$slr6_vs_slr8, model.pair$slr6_vs_slr9, model.pair$slr7_vs_slr8, model.pair$slr7_vs_slr9, model.pair$slr8_vs_slr9, .id = "Comparison")
model.pair.out
write.csv(model.pair.out, file = "PERMANOVA_area_model_outputs.csv")

#------------------------------------
####Kruskal-Wallis####

# Kruskal-Wallis tests for each shape
square.kw <- compare_means(square.sqrt ~ location, data=prism, method = "kruskal.test")
square.kw
# A tibble: 1 x 6
# .y.        p p.adj p.format p.signif method        
# <chr>  <dbl> <dbl> <chr>    <chr>    <chr>         
#   1 square 0.102   0.1 0.1      ns       Kruskal-Wallis

rectangle.kw <- compare_means(rectangle.sqrt ~ location, data=prism, method = "kruskal.test")
rectangle.kw
# # A tibble: 1 x 6
# .y.               p    p.adj p.format p.signif method        
# <chr>         <dbl>    <dbl> <chr>    <chr>    <chr>         
#   1 rectangle 0.0000442 0.000044 4.4e-05  ****     Kruskal-Wallis

circle.kw <- compare_means(circle.sqrt ~ location, data=prism, method = "kruskal.test")
circle.kw
# # A tibble: 1 x 6
# .y.        p p.adj p.format p.signif method        
# <chr>  <dbl> <dbl> <chr>    <chr>    <chr>         
#   1 circle 0.107  0.11 0.11     ns       Kruskal-Wallis

kruskal <- rbind(square.kw,rectangle.kw,circle.kw)
write.csv(kruskal, file="Kruskal_outputs.csv")

# pairwise Dunn tests by bank
square.dunn <- dunnTest(square.sqrt ~ location, data=prism, method = "bh")$res
rectangle.dunn <- dunnTest(rectangle.sqrt ~ location, data=prism, method = "bh")$res
circle.dunn <- dunnTest(circle.sqrt ~ location, data=prism, method = "bh")$res

dunn <- rbind(square.dunn,rectangle.dunn,circle.dunn)
write.csv(dunn, file="Dunn_outputs.csv")

# removing template comparisons in rectangle.dunn for plotting
# not needed for square or circle since no sig diff
rectangle.dunn <- rectangle.dunn[-c(4:6),]
rectangle.dunn

# compact letter display for figures
square.list <- cldList(P.adj ~ Comparison, data = square.dunn, threshold = 0.05)
rectangle.list <- cldList(P.adj ~ Comparison, data = rectangle.dunn, threshold = 0.05)
circle.list <- cldList(P.adj ~ Comparison, data = circle.dunn, threshold = 0.05)

#-----------------------------------
####plots####

# subsetting dataset to remove template
# will instead use line at template value in plots
prism.sub <- prism[ which(prism$location!='Template'), ]
head(prism.sub)

prism.sub$location=factor(prism.sub$location, levels=unique(prism.sub$location)) 

# boxplots comparing shape areas to template
square.box <-
  ggboxplot(
    prism.sub,
    x = "location",
    y = "square",
    color = "grey30",
    palette = c("grey", "#3070cf", "#79d7fb"),
    fill = "location",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Location",
           y = "Area (cm2)",
           title = "Square",
           fill = 'Location') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=2,hjust=0.5, vjust=4) + 
  geom_hline(yintercept = 40.3225, color = "black", size=1) 
# geom_text(data=square.list, aes(x = Group, y=0, vjust=-2.5, label=Letter))
square.box

rectangle.box <-
  ggboxplot(
    prism.sub,
    x = "location",
    y = "rectangle",
    color = "grey30",
    palette = c("grey", "#3070cf", "#79d7fb"),
    fill = "location",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Location",
           y = element_blank(),
           title = "Rectangle",
           fill = 'Location') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=2,hjust=0.5, vjust=4) + 
  geom_hline(yintercept = 12.9032, color = "black", size=1) +
  geom_text(data=rectangle.list, aes(x = Group, y=23, label=Letter)) +
  ylim(8,26)
rectangle.box

circle.box <-
  ggboxplot(
    prism.sub,
    x = "location",
    y = "circle",
    color = "grey30",
    palette = c("grey", "#3070cf", "#79d7fb"),
    fill = "location",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Location",
           y = element_blank(),
           title = "Circle",
           fill = 'Location') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=2,hjust=0.5, vjust=4) + 
  geom_hline(yintercept = 45.60367312, color = "black", size=1) 
# geom_text(data=circle.list, aes(x = Group, y=0, vjust=-2.5, label=Letter))
circle.box

plot=grid.arrange(square.box, rectangle.box, circle.box, ncol=3, nrow=1, widths=c(3.1,3,3), heights=c(3))

#saves plot as PDF
ggsave("prism_area_means.pdf", plot= plot, width=12, height=4, units="in", dpi=300)
# following export, need to change the p values in these plots to reflect the original stats
# the exported stats are based on the subset dataframe, not including template
# the earlier stats include the template in comparisons

#-------------------------------
####shape error normality####

# repeating the same stats and plots with the error measurements

shapiro.test(prism.sub$square.diff)
# Shapiro-Wilk normality test
# 
# data:  prism.sub$square.diff
# W = 0.76877, p-value < 2.2e-16

shapiro.test(prism.sub$rectangle.diff)
# Shapiro-Wilk normality test
# 
# data:  prism.sub$rectangle.diff
# W = 0.66618, p-value < 2.2e-16

shapiro.test(prism.sub$circle.diff)
# Shapiro-Wilk normality test
# 
# data:  prism.sub$circle.diff
# W = 0.77602, p-value = 2.904e-09

# very non-normal, appyling sqrt transformation

prism.sub$square.diff.sqrt <- sqrt(prism.sub$square.diff)
prism.sub$rectangle.diff.sqrt <- sqrt(prism.sub$rectangle.diff)
prism.sub$circle.diff.sqrt <- sqrt(prism.sub$circle.diff)
head(prism.sub)

shapiro.test(prism.sub$square.diff.sqrt)
# Shapiro-Wilk normality test
# 
# data:  prism.sub$square.diff.sqrt
# W = 0.95065, p-value = 1.908e-08

shapiro.test(prism.sub$rectangle.diff.sqrt)
# Shapiro-Wilk normality test
# 
# data:  prism.sub$rectangle.diff.sqrt
# W = 0.9146, p-value = 7.662e-12

shapiro.test(prism.sub$circle.diff.sqrt)
# Shapiro-Wilk normality test
# 
# data:  prism.sub$circle.diff.sqrt
# W = 0.92993, p-value = 0.000506

# still not normal

#------------------------------------
####shape error PERMANOVA####

# Setting seed allows randomized processes to be repeated later
set.seed(215)

# creates a dissimilarity matrix of square-root transformed shape areas, based on Euclidean distance
error.dist <- vegdist(prism.sub[c(17:19)], method="euclidean", na.rm = TRUE)
# replaces missing data (NAs) with 0
error.dist[is.na(error.dist)] <- 0

# testing for homogeneity of variance among locations
error.disp <- betadisper(error.dist, group=prism.sub$location)
permutest(error.disp, bias.adjust = TRUE, perm = 9999)
# significant test value indicates heterogeneous variance

# post hoc tests
error.disp.HSD <- TukeyHSD(error.disp)
error.disp.HSD
# indicates unequal variance in site SLR

# plots showing higher variance in SLR versus other locations
plot(error.disp, ellipse = TRUE, hull = FALSE)
boxplot(error.disp)
# These results suggest unequal variance among locations, with highest variance attributed to SLR
# Since SLR is also the smallest group, this can result in liberal PERMANOVA results (per Anderson and Walsh 2013)
# As a result, the highly significant model terms in the PERMANOVA below should be interpreted with caution

# running the PERMANOVA
# model is nested within location using '/', and is accounted as a random factor using 'strata'
error.perm <- adonis(error.dist ~ location*location/model*face, data=prism.sub, strata = prism.sub$model, permutations = 1e6, parallel = getOption("mc.cores"))
error.perm

write.csv(error.perm$aov.tab, file = "PERMANOVA_error_outputs.csv")
# Significant differences among all main factors

# pairwise PERMANOVA based on significant factors
location.error.pair <- pairwise.adonis2(error.dist ~ location, data = prism.sub, strata = 'model', nperm = 1e6)
location.error.pair

location.error.pair.out <- bind_rows(location.error.pair$Pool_vs_LBTS, location.error.pair$Pool_vs_SLR, location.error.pair$LBTS_vs_SLR, .id = "Comparison")
location.error.pair.out
write.csv(location.error.pair.out, file = "PERMANOVA_error_location_outputs.csv")

model.error.pair <- pairwise.adonis2(error.dist ~ model, data = prism.sub, strata = 'model', nperm = 1e6)
model.error.pair

model.error.pair.out <- bind_rows(model.error.pair$pool1_vs_pool3, model.error.pair$pool1_vs_pool4, model.error.pair$pool1_vs_pool5, model.error.pair$pool1_vs_pool6, model.error.pair$pool1_vs_pool7, model.error.pair$pool3_vs_pool4, model.error.pair$pool3_vs_pool5, model.error.pair$pool3_vs_pool6, model.error.pair$pool3_vs_pool7, model.error.pair$pool4_vs_pool5, model.error.pair$pool4_vs_pool6, model.error.pair$pool4_vs_pool7, model.error.pair$pool5_vs_pool6, model.error.pair$pool5_vs_pool7, model.error.pair$pool6_vs_pool7, model.error.pair$lbts1_vs_lbts2, model.error.pair$lbts1_vs_lbts3, model.error.pair$lbts1_vs_lbts5, model.error.pair$lbts1_vs_lbts6, model.error.pair$lbts1_vs_lbts7, model.error.pair$lbts1_vs_lbts8, model.error.pair$lbts2_vs_lbts3, model.error.pair$lbts2_vs_lbts5, model.error.pair$lbts2_vs_lbts6, model.error.pair$lbts2_vs_lbts7, model.error.pair$lbts2_vs_lbts8, model.error.pair$lbts3_vs_lbts5, model.error.pair$lbts3_vs_lbts6, model.error.pair$lbts3_vs_lbts7, model.error.pair$lbts3_vs_lbts8, model.error.pair$lbts5_vs_lbts6, model.error.pair$lbts5_vs_lbts7, model.error.pair$lbts5_vs_lbts8, model.error.pair$lbts6_vs_lbts7, model.error.pair$lbts6_vs_lbts8, model.error.pair$lbts7_vs_lbts8, model.error.pair$slr1_vs_slr2, model.error.pair$slr1_vs_slr6, model.error.pair$slr1_vs_slr7, model.error.pair$slr1_vs_slr8, model.error.pair$slr1_vs_slr9, model.error.pair$slr2_vs_slr6, model.error.pair$slr2_vs_slr7, model.error.pair$slr2_vs_slr8, model.error.pair$slr2_vs_slr9, model.error.pair$slr6_vs_slr7, model.error.pair$slr6_vs_slr8, model.error.pair$slr6_vs_slr9, model.error.pair$slr7_vs_slr8, model.error.pair$slr7_vs_slr9, model.error.pair$slr8_vs_slr9, .id = "Comparison")
model.error.pair.out
write.csv(model.error.pair.out, file = "PERMANOVA_error_model_outputs.csv")

face.error.pair <- pairwise.adonis2(error.dist ~ face, data = prism.sub, strata = 'model', nperm = 1e6)
face.error.pair

face.error.pair.out <- bind_rows(face.error.pair$side1_vs_side2, face.error.pair$side1_vs_side3, face.error.pair$side1_vs_side4, face.error.pair$side1_vs_top, face.error.pair$side2_vs_side3, face.error.pair$side2_vs_side4, face.error.pair$side2_vs_top, face.error.pair$side3_vs_side4, face.error.pair$side3_vs_top, face.error.pair$side4_vs_top, .id = "Comparison")
face.error.pair.out
write.csv(face.error.pair.out, file = "PERMANOVA_error_face_outputs.csv")

#------------------------------------
####exporting/importing results####

# saving dataframes and test results
save(model, model.anova, prism, area.dist, area.perm, model.pair, prism.sub, error.dist, error.perm, location.error.pair, model.error.pair, face.error.pair, file = "prism_stats.RData")
# to load back in
load("prism_stats.RData")

#------------------------------------
####shape error Kruskal-Wallis####

# Kruskal-Wallis tests for each shape
square.diff.kw <- compare_means(square.diff.sqrt ~ location, data=prism.sub, method = "kruskal.test")
square.diff.kw
# # A tibble: 1 x 6
# .y.                        p      p.adj p.format p.signif method        
# <chr>                  <dbl>      <dbl> <chr>    <chr>    <chr>         
#   1 square.diff.sqrt 0.000000311 0.00000031 3.1e-07  ****     Kruskal-Wallis

rectangle.diff.kw <- compare_means(rectangle.diff.sqrt ~ location, data=prism.sub, method = "kruskal.test")
rectangle.diff.kw
# # A tibble: 1 x 6
# .y.                          p     p.adj p.format p.signif method        
# <chr>                    <dbl>     <dbl> <chr>    <chr>    <chr>         
#   1 rectangle.diff.sqrt 0.00000453 0.0000045 4.5e-06  ****     Kruskal-Wallis

circle.diff.kw <- compare_means(circle.diff.sqrt ~ location, data=prism.sub, method = "kruskal.test")
circle.diff.kw
# # A tibble: 1 x 6
# .y.                    p  p.adj p.format p.signif method        
# <chr>              <dbl>  <dbl> <chr>    <chr>    <chr>         
#   1 circle.diff.sqrt 0.00631 0.0063 0.0063   **       Kruskal-Wallis

kruskal <- rbind(square.diff.kw,rectangle.diff.kw,circle.diff.kw)
write.csv(kruskal, file="Kruskal_error_outputs.csv")

# pairwise Dunn tests by location
square.diff.dunn <- dunnTest(square.diff.sqrt ~ location, data=prism.sub, method = "bh")$res
rectangle.diff.dunn <- dunnTest(rectangle.diff.sqrt ~ location, data=prism.sub, method = "bh")$res
circle.diff.dunn <- dunnTest(circle.diff.sqrt ~ location, data=prism.sub, method = "bh")$res

dunn <- rbind(square.diff.dunn,rectangle.diff.dunn,circle.diff.dunn)
write.csv(dunn, file="Dunn_error_outputs.csv")

# compact letter display for figures
square.diff.list <- cldList(P.adj ~ Comparison, data = square.diff.dunn, threshold = 0.05)
rectangle.diff.list <- cldList(P.adj ~ Comparison, data = rectangle.diff.dunn, threshold = 0.05)
circle.diff.list <- cldList(P.adj ~ Comparison, data = circle.diff.dunn, threshold = 0.05)

#-----------------------------------
####shape error plots####

prism.sub$location=factor(prism.sub$location, levels=unique(prism.sub$location)) 

# creates a function that forces y axis labels to have 1 decimal place
scale <- function(x) sprintf("%.0f", x)

# boxplots comparing shape areas to template
square.diff.box <-
  ggboxplot(
    prism.sub,
    x = "location",
    y = "square.diff",
    color = "grey30",
    palette = c("grey", "#3070cf", "#79d7fb"),
    fill = "location",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Location",
           y = "Shape Error (cm2)",
           #title = "Square",
           fill = 'Location') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none") + 
  stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=2,hjust=0.5, vjust=4) +
 geom_text(data=square.diff.list, aes(x = Group, y=0, vjust=-22.5, label=Letter)) +
 scale_y_continuous(labels = scale)
square.diff.box

rectangle.diff.box <-
  ggboxplot(
    prism.sub,
    x = "location",
    y = "rectangle.diff",
    color = "grey30",
    palette = c("grey", "#3070cf", "#79d7fb"),
    fill = "location",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Location",
           y = element_blank(),
           #title = "Rectangle",
           fill = 'Location') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none") + 
  stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=2,hjust=0.5, vjust=4) +
  geom_text(data=rectangle.diff.list, aes(x = Group, y=0, vjust=-22.5, label=Letter))  +
  scale_y_continuous(labels = scale)
rectangle.diff.box

circle.diff.box <-
  ggboxplot(
    prism.sub,
    x = "location",
    y = "circle.diff",
    color = "grey30",
    palette = c("grey", "#3070cf", "#79d7fb"),
    fill = "location",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Location",
           y = element_blank(),
           #title = "Circle",
           fill = 'Location') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none") + 
  stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=2,hjust=0.5, vjust=4) + 
  geom_text(data=circle.diff.list, aes(x = Group, y=0, vjust=-22.5, label=Letter)) +
  scale_y_continuous(labels = scale)
circle.diff.box

plot2=grid.arrange(square.diff.box, rectangle.diff.box, circle.diff.box, ncol=3, nrow=1, widths=c(3.1,3,3), heights=c(3))

#saves plot as PDF
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
