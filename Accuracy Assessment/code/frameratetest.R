

frame <- read.csv("../data/frameratetest.csv", header = TRUE)


head(frame)



##################################################################
 
#Anova on frame rate

frame.kw <- kruskal.test(total ~ frame.rate, data = frame)
frame.kw


