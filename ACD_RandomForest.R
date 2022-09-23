### ANT CHCS AND DESICCATION - RANDOM FOREST SELECTION OF CHCS ###
## BLOCK 0: Set up + clean up -------------------------------------------------------

## Housekeeping
cat("\014") # Clear console
rm(list=ls()) # Remove all variables
options(scipen = 1, digits=4) # scientific notation setting

## Libraries
library(ggplot2) # plotting
library(ggpubr) # for ggscatter plots
library(GGally) # for ggpairs correlation plots
library(caTools) # for sampling the data set into training and test sets
library(randomForest) # for the random forest algorithm
library(dplyr) # data transforming

## Set WD + load data
#setwd("F:/Users/Brian/Google Drive (ba.whyte@berkeley.edu)/R_folder/1 - NeilNSF/1 - Publication")
setwd("H:/My Drive/0 - R/Ant-CHCs-Desiccation")
#setwd("/Volumes/GoogleDrive/My Drive/0 - R/1 - NeilNSF/1 - Publication/Resubmitting")
data <- read.csv("ModelData_byTube_Sep2022.csv")

df1 <- subset(data, select = c(17:88)) # only individual CHCs
bm <- data$Avg
df2 <- df1 / t(bm)
df <- cbind(data$LT50d, df2)
names(df)[names(df) == 'data$LT50d'] <- 'LT50'
#

## BLOCK 1: Example using Iris dataset -------------------------------------------------------
#https://www.geeksforgeeks.org/random-forest-approach-in-r-programming/
# Load Iris
data(iris) # loads the data set which can be called by this same name
# Splitting data in train and test data
split <- sample.split(iris, SplitRatio = 0.7) # randomly chooses which columns to select?
train <- subset(iris, split == "TRUE") # 70% of iris data
test <- subset(iris, split == "FALSE") # 30% of iris data

set.seed(120)  # Setting seed. IIRC this is preferred before randomization in general
classifier_RF = randomForest(x = train[-5], # the -5 removes the species character column
                             y = train$Species, # because species will be the Y we try to predict?
                             ntree = 500) # num. of trees to generate
# Predicting the Test set results
y_pred = predict(classifier_RF, newdata = test[-5]) # also remove species column from test dataset
# Confusion Matrix
confusion_mtx = table(test[, 5], y_pred) # This is similar to the confusion matrix output in classifier_RF
# Plotting model
plot(classifier_RF) # I guess this shows how the error levels out after enough trees generated?
# Importance plot
importance(classifier_RF) # For the raw importance inference
# Variable importance plot
varImpPlot(classifier_RF) # THIS ONE is great. A better visualization of the same importance result. 

## BLOCK 2: Example but with individual CHCs to predict LT50d -------------------------------------------------------
# Split the data set into training and test sets
split <- sample.split(df, SplitRatio = 0.7)
train <- subset(df, split == "TRUE") # 70% of data
test <- subset(df, split == "FALSE") # 30% of data

set.seed(120)  # Setting seed. IIRC this is preferred before randomization in general
reg_RF = randomForest(x = train[-1], # remove the LT50d column
                             y = train$LT50, # because LT50d will be the Y we try to predict
                             ntree = 1000, type = "regression")
y_pred = predict(reg_RF, newdata = test[-1]) # also remove species column from test dataset
confusion_mtx = table(test[,1], y_pred)
#importance(reg_RF) # showing only top 30 variables, I think
#varImpPlot(reg_RF)
v <- as.data.frame(importance(reg_RF))
v$new <- rownames(v)
top <- v[order(-v$IncNodePurity),]
head(top,10) # view top 10 on ordered list
#write.csv(top5, file = "IncNodePurityList.csv")
t <- c() # empty list to hold all top 10 CHC names
for (i in 1:100) {
  reg_RF = randomForest(x = train[-1], # remove the LT50d column
                        y = train$LT50, # because LT50d will be the Y we try to predict
                        ntree = 1000, type = "regression")
  v <- as.data.frame(importance(reg_RF))
  v$new <- rownames(v)
  top <- v[order(-v$IncNodePurity),]
  t <- c(t,top$new[1:10])
  #print("Done!")
}
plot(as.factor(t))


#varImpPlot(reg_RF)
#varImpPlot(classifier_RF, n.var=(nrow(df))) # showing all variables

## BLOCK 3: Correlating individual CHCs vs. LT50 -------------------------------------------------------------------
# Spearman correlation X vs. Y
g1 <- ggscatter(data, x = "SumCHC", y = "Avg",
                 add = "reg.line", #conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "spearman")

# Correlated pairs plot
#dfp <- subset(df, select = c("LT50","MeC34.3492","n.C33.3300","MeC32.3258","MeC33.3345","TriMeC31.3225"))
dfp <- subset(df, select = c("LT50",rownames(top5)[1:5]))
ggpairs(dfp, lower = list(continuous = wrap("smooth")), upper = list(continous = wrap("cor", method = "spearman")))

# P-values of 72 Spearman correlations between LT50 and each CHC
#df <- cbind(data$LT50d, df2)
#ames(df)[names(df) == 'data$LT50d'] <- 'LT50'
plist <- c() # empty list to store p values from spearman tests
clist <- c() # empty list to store correlation coefficients from spearman tests
for (i in 1:72) {
  y <- df$LT50 # same Y variable each time
  x <- df[,i+1] # current CHC column
  p <- as.numeric(cor.test(x,y,method="spearman",exact=FALSE)[3]) # p-value of current spearman test
  c <- as.numeric(cor.test(x,y,method="spearman",exact=FALSE)[4]) # p-value of current spearman test
  plist <- c(plist,p)
  clist <- c(clist,c)
}

# Spearman data frame using df1 CHC values (CHC mass)
sdf1 <- tibble(colnames(df1),p.adjust(plist,method="BH"),clist)
colnames(sdf1) <- c("CHC","p-adjusted","coef")
sdf1 <- sdf1[order(sdf1$`p-adjusted`),] # lowest value is most significant

# Spearman data frame using df2 CHC values (CHC mass per body mass)
sdf2 <- tibble(colnames(df2),p.adjust(plist,method="BH"),clist)
colnames(sdf2) <- c("CHC","p-adjusted","coef")
sdf2 <- sdf2[order(sdf2$`p-adjusted`),] # lowest value is most significant




