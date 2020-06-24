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
# Df Sum Sq Mean Sq F value   Pr(>F)    
# location     2  65.66   32.83   18.21 7.54e-05 ***
#   Residuals   16  28.85    1.80                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

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

# transfer from Primer to R in progress

# replacing missing data (NA) with '0' for dissimilarity matrix
prism.na <- prism
prism.na[is.na(prism.na)] <- 0

# Setting seed allows randomized processes to be repeated later
set.seed(215)

# creates a dissimilarity matrix of square-root transformed shape areas, based on Euclidean distance
area.dist <- vegdist(prism.na[c(14:16)], method="euclidean")

# running the PERMANOVA
area.perm <- adonis(area.dist ~ location*location/model*face, data=prism.na, permutations = 999, parallel = getOption("mc.cores"))
area.perm
# Call:
#   adonis(formula = area.dist ~ location * location/model * face,      data = prism.na, permutations = 999, parallel = getOption("mc.cores")) 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# location              3     497.3 165.753 12.1457 0.07601  0.001 ***
#   face                  4     285.2  71.290  5.2239 0.04359  0.001 ***
#   location:model       16     258.0  16.126  1.1816 0.03944  0.226    
# location:face         8     178.6  22.324  1.6358 0.02730  0.052 .  
# location:model:face  64    1433.8  22.403  1.6416 0.21916  0.001 ***
#   Residuals           285    3889.4  13.647         0.59451           
# Total               380    6542.2                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Significant differences between location and prism face, but not among models within location

# pairwise PERMANOVA based on significant factors
loc.perm <- pairwise.adonis(prism.na[c(14:16)], factors = prism.na$location, sim.function = 'vegdist', sim.method = 'euclidean', p.adjust.m = 'bonferroni')
loc.perm
# pairs total.DF    F.Model          R2 p.value p.adjusted sig
# 1 Template vs Pool      120  2.7429921 0.022531006   0.158      0.948    
# 2 Template vs LBTS      140  1.7564323 0.012478523   0.174      1.000    
# 3  Template vs SLR      120  2.0002233 0.016530741   0.137      0.822    
# 4     Pool vs LBTS      259 21.9030524 0.078252281   0.001      0.006   *
#   5      Pool vs SLR      239 26.5216988 0.100262848   0.001      0.006   *
#   6      LBTS vs SLR      259  0.7251707 0.002802861   0.488      1.000    

face.perm <- pairwise.adonis(prism.na[c(14:16)], factors = prism.na$face, sim.function = 'vegdist', sim.method = 'euclidean', p.adjust.m = 'bonferroni')
face.perm
# pairs total.DF     F.Model           R2 p.value p.adjusted sig
# 1    template vs top       76  2.66628501 0.0343300187   0.211      1.000    
# 2  template vs side1       76  1.80902467 0.0235522411   0.319      1.000    
# 3  template vs side2       76  1.77630781 0.0231361453   0.150      1.000    
# 4  template vs side3       76  1.91050527 0.0248406282   0.142      1.000    
# 5  template vs side4       76  1.80571931 0.0235102193   0.172      1.000    
# 6       top vs side1      151  9.27424166 0.0582281326   0.002      0.030   .
# 7       top vs side2      151 13.58750611 0.0830595590   0.001      0.015   .
# 8       top vs side3      151 14.35470130 0.0873397669   0.001      0.015   .
# 9       top vs side4      151 11.46005165 0.0709776290   0.001      0.015   .
# 10    side1 vs side2      151  0.40782679 0.0027114732   0.693      1.000    
# 11    side1 vs side3      151  0.46885095 0.0031159336   0.643      1.000    
# 12    side1 vs side4      151  0.33394118 0.0022213293   0.759      1.000    
# 13    side2 vs side3      151  0.02519249 0.0001679217   0.989      1.000    
# 14    side2 vs side4      151  0.35821957 0.0023824409   0.722      1.000    
# 15    side3 vs side4      151  0.39017569 0.0025944227   0.709      1.000    

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

# transfer from Primer to R in progress

# replacing missing data (NA) with '0' for dissimilarity matrix
prism.sub.na <- prism.sub
prism.sub.na[is.na(prism.sub.na)] <- 0

# Setting seed allows randomized processes to be repeated later
set.seed(215)

# creates a dissimilarity matrix of square-root transformed shape areas, based on Euclidean distance
error.dist <- vegdist(prism.sub.na[c(17:19)], method="euclidean")

# running the PERMANOVA
error.perm <- adonis(error.dist ~ location*location/model*face, data=prism.sub.na, permutations = 999, parallel = getOption("mc.cores"))
error.perm
# Call:
#   adonis(formula = error.dist ~ location * location/model * face,      data = prism.sub.na, permutations = 999, parallel = getOption("mc.cores")) 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# location              2      9.09  4.5455  3.4104 0.01309  0.007 ** 
#   face                  4     10.84  2.7109  2.0339 0.01561  0.038 *  
#   location:model       16     75.18  4.6986  3.5252 0.10821  0.001 ***
#   location:face         8     16.77  2.0962  1.5727 0.02414  0.060 .  
# location:model:face  64    203.02  3.1722  2.3800 0.29222  0.001 ***
#   Residuals           285    379.86  1.3329         0.54675           
# Total               379    694.77                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Significant differences among all main factors

# pairwise PERMANOVA based on significant factors
loc.error.perm <- pairwise.adonis(prism.sub.na[c(17:19)], factors = prism.sub.na$location, sim.function = 'vegdist', sim.method = 'euclidean', p.adjust.m = 'bonferroni')
loc.error.perm
# pairs total.DF   F.Model           R2 p.value p.adjusted sig
# 1 Pool vs LBTS      259 3.8379523 0.0146577388   0.014      0.042   .
# 2  Pool vs SLR      239 0.2094955 0.0008794589   0.870      1.000    
# 3  LBTS vs SLR      259 3.3704834 0.0128954246   0.037      0.111    

face.error.perm <- pairwise.adonis(prism.sub.na[c(17:19)], factors = prism.sub.na$face, sim.function = 'vegdist', sim.method = 'euclidean', p.adjust.m = 'bonferroni')
face.error.perm
# pairs total.DF   F.Model           R2 p.value p.adjusted sig
# 1    top vs side1      151 0.2291161 0.0015251110   0.854       1.00    
# 2    top vs side2      151 4.0068074 0.0260170800   0.012       0.12    
# 3    top vs side3      151 3.7895076 0.0246408723   0.020       0.20    
# 4    top vs side4      151 1.8375157 0.0121018558   0.155       1.00    
# 5  side1 vs side2      151 2.1879826 0.0143768418   0.099       0.99    
# 6  side1 vs side3      151 1.9100300 0.0125734290   0.126       1.00    
# 7  side1 vs side4      151 0.8543645 0.0056635055   0.448       1.00    
# 8  side2 vs side3      151 0.0861483 0.0005739924   0.958       1.00    
# 9  side2 vs side4      151 0.4092606 0.0027209801   0.715       1.00    
# 10 side3 vs side4      151 0.2461135 0.0016380690   0.836       1.00   

model.error.perm <- pairwise.adonis(prism.sub.na[c(17:19)], factors = prism.sub.na$model, sim.function = 'vegdist', sim.method = 'euclidean', p.adjust.m = 'bonferroni')
model.error.perm
# pairs total.DF      F.Model           R2 p.value p.adjusted sig
# 1   pool1 vs pool3       39  4.066942928 9.667788e-02   0.009      1.000    
# 2   pool1 vs pool4       39  1.701360186 4.285395e-02   0.178      1.000    
# 3   pool1 vs pool5       39  5.804774556 1.325147e-01   0.003      0.513    
# 4   pool1 vs pool6       39 10.146247519 2.107381e-01   0.001      0.171    
# 5   pool1 vs pool7       39  2.524522597 6.229617e-02   0.073      1.000    
# 6   pool1 vs lbts1       39  1.363324122 3.463437e-02   0.267      1.000    
# 7   pool1 vs lbts2       39  2.058600484 5.138973e-02   0.111      1.000    
# 8   pool1 vs lbts3       39  3.421000587 8.259097e-02   0.019      1.000    
# 9   pool1 vs lbts5       39  5.015379080 1.165950e-01   0.002      0.342    
# 10  pool1 vs lbts6       39 10.889323896 2.227342e-01   0.001      0.171    
# 11  pool1 vs lbts7       39  1.094227398 2.798949e-02   0.307      1.000    
# 12  pool1 vs lbts8       39  1.102911428 2.820535e-02   0.371      1.000    
# 13   pool1 vs slr1       39  4.474616817 1.053480e-01   0.010      1.000    
# 14   pool1 vs slr2       39  1.875267371 4.702833e-02   0.169      1.000    
# 15   pool1 vs slr6       39  4.605924790 1.081053e-01   0.011      1.000    
# 16   pool1 vs slr7       39  1.585787798 4.005952e-02   0.201      1.000    
# 17   pool1 vs slr8       39  3.121713319 7.591399e-02   0.039      1.000    
# 18   pool1 vs slr9       39  0.624653456 1.617240e-02   0.586      1.000    
# 19  pool3 vs pool4       39  0.540838451 1.403287e-02   0.631      1.000    
# 20  pool3 vs pool5       39  8.477204226 1.823949e-01   0.001      0.171    
# 21  pool3 vs pool6       39 12.475029668 2.471525e-01   0.001      0.171    
# 22  pool3 vs pool7       39  6.440979866 1.449333e-01   0.003      0.513    
# 23  pool3 vs lbts1       39  1.852594236 4.648616e-02   0.155      1.000    
# 24  pool3 vs lbts2       39  4.990169635 1.160770e-01   0.005      0.855    
# 25  pool3 vs lbts3       39  3.010029034 7.339739e-02   0.053      1.000    
# 26  pool3 vs lbts5       39  7.611252107 1.668722e-01   0.001      0.171    
# 27  pool3 vs lbts6       39 14.774016007 2.799487e-01   0.001      0.171    
# 28  pool3 vs lbts7       39  5.078424405 1.178879e-01   0.005      0.855    
# 29  pool3 vs lbts8       39  2.614721665 6.437867e-02   0.054      1.000    
# 30   pool3 vs slr1       39  8.801288950 1.880566e-01   0.001      0.171    
# 31   pool3 vs slr2       39  6.732637071 1.505084e-01   0.003      0.513    
# 32   pool3 vs slr6       39 10.493257724 2.163859e-01   0.001      0.171    
# 33   pool3 vs slr7       39  0.003545494 9.329378e-05   1.000      1.000    
# 34   pool3 vs slr8       39  0.183848167 4.814815e-03   0.868      1.000    
# 35   pool3 vs slr9       39  1.858300902 4.662268e-02   0.128      1.000    
# 36  pool4 vs pool5       39  5.591053399 1.282615e-01   0.002      0.342    
# 37  pool4 vs pool6       39  9.121451837 1.935732e-01   0.001      0.171    
# 38  pool4 vs pool7       39  3.545579611 8.534192e-02   0.019      1.000    
# 39  pool4 vs lbts1       39  0.723232335 1.867696e-02   0.531      1.000    
# 40  pool4 vs lbts2       39  2.788306038 6.836043e-02   0.052      1.000    
# 41  pool4 vs lbts3       39  2.069173497 5.164003e-02   0.124      1.000    
# 42  pool4 vs lbts5       39  5.220905632 1.207958e-01   0.004      0.684    
# 43  pool4 vs lbts6       39 11.286107041 2.289917e-01   0.001      0.171    
# 44  pool4 vs lbts7       39  2.736718627 6.718064e-02   0.070      1.000    
# 45  pool4 vs lbts8       39  1.263771862 3.218672e-02   0.310      1.000    
# 46   pool4 vs slr1       39  5.955634313 1.354919e-01   0.006      1.000    
# 47   pool4 vs slr2       39  4.010986394 9.547470e-02   0.019      1.000    
# 48   pool4 vs slr6       39  7.231757682 1.598823e-01   0.001      0.171    
# 49   pool4 vs slr7       39  0.229181348 5.994932e-03   0.827      1.000    
# 50   pool4 vs slr8       39  0.829165179 2.135419e-02   0.459      1.000    
# 51   pool4 vs slr9       39  0.485827813 1.262355e-02   0.698      1.000    
# 52  pool5 vs pool6       39  0.541401250 1.404726e-02   0.649      1.000    
# 53  pool5 vs pool7       39  0.802035164 2.066993e-02   0.467      1.000    
# 54  pool5 vs lbts1       39  1.459362361 3.698393e-02   0.237      1.000    
# 55  pool5 vs lbts2       39  0.367775115 9.585521e-03   0.799      1.000    
# 56  pool5 vs lbts3       39  1.272022896 3.239005e-02   0.290      1.000    
# 57  pool5 vs lbts5       39  0.321644606 8.393288e-03   0.801      1.000    
# 58  pool5 vs lbts6       39  2.176662606 5.417729e-02   0.113      1.000    
# 59  pool5 vs lbts7       39  0.789589832 2.035571e-02   0.461      1.000    
# 60  pool5 vs lbts8       39  0.426512117 1.109942e-02   0.740      1.000    
# 61   pool5 vs slr1       39  0.472840830 1.229025e-02   0.670      1.000    
# 62   pool5 vs slr2       39  0.928279698 2.384590e-02   0.390      1.000    
# 63   pool5 vs slr6       39  2.143240333 5.338982e-02   0.084      1.000    
# 64   pool5 vs slr7       39  3.581976224 8.614252e-02   0.039      1.000    
# 65   pool5 vs slr8       39  6.827570808 1.523074e-01   0.002      0.342    
# 66   pool5 vs slr9       39  2.662306396 6.547357e-02   0.076      1.000    
# 67  pool6 vs pool7       39  2.497854754 6.167869e-02   0.055      1.000    
# 68  pool6 vs lbts1       39  3.041486716 7.410762e-02   0.025      1.000    
# 69  pool6 vs lbts2       39  1.155105348 2.950076e-02   0.317      1.000    
# 70  pool6 vs lbts3       39  2.406920205 5.956703e-02   0.063      1.000    
# 71  pool6 vs lbts5       39  0.449527346 1.169136e-02   0.689      1.000    
# 72  pool6 vs lbts6       39  0.982357692 2.520006e-02   0.389      1.000    
# 73  pool6 vs lbts7       39  1.348667702 3.427480e-02   0.234      1.000    
# 74  pool6 vs lbts8       39  1.173602470 2.995901e-02   0.347      1.000    
# 75   pool6 vs slr1       39  0.107208794 2.813347e-03   0.943      1.000    
# 76   pool6 vs slr2       39  1.147186378 2.930444e-02   0.307      1.000    
# 77   pool6 vs slr6       39  1.733006536 4.361630e-02   0.159      1.000    
# 78   pool6 vs slr7       39  5.019533909 1.166803e-01   0.015      1.000    
# 79   pool6 vs slr8       39  9.222830285 1.953045e-01   0.001      0.171    
# 80   pool6 vs slr9       39  4.913950489 1.145071e-01   0.012      1.000    
# 81  pool7 vs lbts1       39  0.835524702 2.151444e-02   0.490      1.000    
# 82  pool7 vs lbts2       39  0.199129946 5.212945e-03   0.891      1.000    
# 83  pool7 vs lbts3       39  1.568894118 3.964968e-02   0.216      1.000    
# 84  pool7 vs lbts5       39  1.283498316 3.267271e-02   0.272      1.000    
# 85  pool7 vs lbts6       39  4.351829114 1.027542e-01   0.017      1.000    
# 86  pool7 vs lbts7       39  0.197725531 5.176369e-03   0.846      1.000    
# 87  pool7 vs lbts8       39  0.202090030 5.290025e-03   0.895      1.000    
# 88   pool7 vs slr1       39  1.197390387 3.054771e-02   0.287      1.000    
# 89   pool7 vs slr2       39  0.514920469 1.336938e-02   0.608      1.000    
# 90   pool7 vs slr6       39  2.266471635 5.628682e-02   0.096      1.000    
# 91   pool7 vs slr7       39  2.613319425 6.434636e-02   0.080      1.000    
# 92   pool7 vs slr8       39  5.194916649 1.202669e-01   0.006      1.000    
# 93   pool7 vs slr9       39  1.205168498 3.074004e-02   0.324      1.000    
# 94  lbts1 vs lbts2       39  0.642908394 1.663716e-02   0.582      1.000    
# 95  lbts1 vs lbts3       39  0.467360396 1.214953e-02   0.662      1.000    
# 96  lbts1 vs lbts5       39  1.655207208 4.173997e-02   0.176      1.000    
# 97  lbts1 vs lbts6       39  4.819089045 1.125453e-01   0.008      1.000    
# 98  lbts1 vs lbts7       39  1.076049690 2.753732e-02   0.313      1.000    
# 99  lbts1 vs lbts8       39  0.150367731 3.941449e-03   0.926      1.000    
# 100  lbts1 vs slr1       39  2.340870175 5.802726e-02   0.083      1.000    
# 101  lbts1 vs slr2       39  1.725327558 4.343142e-02   0.164      1.000    
# 102  lbts1 vs slr6       39  3.562359182 8.571119e-02   0.030      1.000    
# 103  lbts1 vs slr7       39  0.963523634 2.472886e-02   0.353      1.000    
# 104  lbts1 vs slr8       39  2.184004853 5.435010e-02   0.105      1.000    
# 105  lbts1 vs slr9       39  0.284244428 7.424580e-03   0.854      1.000    
# 106 lbts2 vs lbts3       39  0.898169993 2.309029e-02   0.442      1.000    
# 107 lbts2 vs lbts5       39  0.493738183 1.282645e-02   0.665      1.000    
# 108 lbts2 vs lbts6       39  2.252843118 5.596730e-02   0.102      1.000    
# 109 lbts2 vs lbts7       39  0.467238150 1.214639e-02   0.637      1.000    
# 110 lbts2 vs lbts8       39  0.102153951 2.681054e-03   0.948      1.000    
# 111  lbts2 vs slr1       39  0.759326461 1.959081e-02   0.497      1.000    
# 112  lbts2 vs slr2       39  0.632996059 1.638486e-02   0.567      1.000    
# 113  lbts2 vs slr6       39  1.479323712 3.747085e-02   0.219      1.000    
# 114  lbts2 vs slr7       39  2.445683947 6.046835e-02   0.078      1.000    
# 115  lbts2 vs slr8       39  4.539721349 1.067172e-01   0.010      1.000    
# 116  lbts2 vs slr9       39  1.229647410 3.134485e-02   0.312      1.000    
# 117 lbts3 vs lbts5       39  1.025100122 2.626771e-02   0.387      1.000    
# 118 lbts3 vs lbts6       39  3.757612248 8.998628e-02   0.016      1.000    
# 119 lbts3 vs lbts7       39  1.918963963 4.807149e-02   0.139      1.000    
# 120 lbts3 vs lbts8       39  0.397710323 1.035766e-02   0.731      1.000    
# 121  lbts3 vs slr1       39  2.152855432 5.361650e-02   0.114      1.000    
# 122  lbts3 vs slr2       39  2.472437951 6.108942e-02   0.078      1.000    
# 123  lbts3 vs slr6       39  3.853985342 9.208168e-02   0.016      1.000    
# 124  lbts3 vs slr7       39  1.720977157 4.332666e-02   0.149      1.000    
# 125  lbts3 vs slr8       39  3.314108522 8.021736e-02   0.036      1.000    
# 126  lbts3 vs slr9       39  1.473039368 3.731761e-02   0.230      1.000    
# 127 lbts5 vs lbts6       39  0.983628987 2.523185e-02   0.400      1.000    
# 128 lbts5 vs lbts7       39  1.317103792 3.349951e-02   0.262      1.000    
# 129 lbts5 vs lbts8       39  0.611556789 1.583870e-02   0.615      1.000    
# 130  lbts5 vs slr1       39  0.416692377 1.084665e-02   0.724      1.000    
# 131  lbts5 vs slr2       39  1.289035245 3.280903e-02   0.322      1.000    
# 132  lbts5 vs slr6       39  1.592899573 4.023195e-02   0.204      1.000    
# 133  lbts5 vs slr7       39  3.820604994 9.135700e-02   0.023      1.000    
# 134  lbts5 vs slr8       39  6.696500323 1.498216e-01   0.002      0.342    
# 135  lbts5 vs slr9       39  3.008356088 7.335959e-02   0.021      1.000    
# 136 lbts6 vs lbts7       39  2.479975777 6.126426e-02   0.100      1.000    
# 137 lbts6 vs lbts8       39  2.276356805 5.651844e-02   0.084      1.000    
# 138  lbts6 vs slr1       39  0.354247016 9.236187e-03   0.755      1.000    
# 139  lbts6 vs slr2       39  1.862105203 4.671367e-02   0.134      1.000    
# 140  lbts6 vs slr6       39  1.054038242 2.698923e-02   0.330      1.000    
# 141  lbts6 vs slr7       39  6.614749094 1.482637e-01   0.009      1.000    
# 142  lbts6 vs slr8       39 11.169297602 2.271600e-01   0.002      0.342 
# 143  lbts6 vs slr9       39  7.0059725 0.155667618   0.002      0.342    
# 144 lbts7 vs lbts8       39  0.4822868 0.012532695   0.668      1.000    
# 145  lbts7 vs slr1       39  0.7946767 0.020484170   0.416      1.000    
# 146  lbts7 vs slr2       39  0.1002824 0.002632066   0.918      1.000    
# 147  lbts7 vs slr6       39  1.1606270 0.029637599   0.288      1.000    
# 148  lbts7 vs slr7       39  2.5068248 0.061886479   0.080      1.000    
# 149  lbts7 vs slr8       39  4.5440585 0.106808298   0.012      1.000    
# 150  lbts7 vs slr9       39  0.9964874 0.025553259   0.356      1.000    
# 151  lbts8 vs slr1       39  0.9928570 0.025462535   0.395      1.000    
# 152  lbts8 vs slr2       39  0.7699795 0.019860198   0.483      1.000    
# 153  lbts8 vs slr6       39  1.7469422 0.043951611   0.167      1.000    
# 154  lbts8 vs slr7       39  1.4855187 0.037621862   0.227      1.000    
# 155  lbts8 vs slr8       39  2.8465212 0.069688216   0.061      1.000    
# 156  lbts8 vs slr9       39  0.4902791 0.012737737   0.689      1.000    
# 157   slr1 vs slr2       39  0.5281593 0.013708396   0.575      1.000    
# 158   slr1 vs slr6       39  0.6034238 0.015631355   0.544      1.000    
# 159   slr1 vs slr7       39  4.3839772 0.103434775   0.027      1.000    
# 160   slr1 vs slr8       39  7.4770411 0.164413536   0.002      0.342    
# 161   slr1 vs slr9       39  3.2348746 0.078449968   0.047      1.000    
# 162   slr2 vs slr6       39  0.6307163 0.016326809   0.524      1.000    
# 163   slr2 vs slr7       39  3.2668468 0.079163956   0.044      1.000    
# 164   slr2 vs slr8       39  5.6919402 0.130274376   0.008      1.000    
# 165   slr2 vs slr9       39  1.7999323 0.045224507   0.150      1.000    
# 166   slr6 vs slr7       39  5.1032474 0.118395891   0.012      1.000    
# 167   slr6 vs slr8       39  8.2643783 0.178633727   0.001      0.171    
# 168   slr6 vs slr9       39  4.2893487 0.101428582   0.019      1.000    
# 169   slr7 vs slr8       39  0.1122000 0.002943940   0.934      1.000    
# 170   slr7 vs slr9       39  0.8631525 0.022210049   0.428      1.000    
# 171   slr8 vs slr9       39  2.0339178 0.050804865   0.120      1.000   

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
