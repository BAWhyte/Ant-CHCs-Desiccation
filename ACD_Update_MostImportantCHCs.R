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
library(pls) # for principle component regression

## Set WD + load data
#setwd("F:/Users/Brian/Google Drive (ba.whyte@berkeley.edu)/R_folder/1 - NeilNSF/1 - Publication")
setwd("H:/My Drive/0 - R/Ant-CHCs-Desiccation")
#setwd("/Volumes/GoogleDrive/My Drive/0 - R/1 - NeilNSF/1 - Publication/Resubmitting")
data <- read.csv("ModelData_byTube_Sep2022.csv")

## Data frame 1 = Raw masses, individual CHCs
df1 <- subset(data, select = c(17:88)) # only individual CHCs
df1 <- cbind(data$LT50d, df1)
names(df1)[names(df2) == 'data$LT50d'] <- 'LT50'
## Data frame 2 = CHC mass per body mass, individual CHCs
df2 <- subset(data, select = c(17:88)) # only individual CHCs
bm <- data$Avg
df2 <- df2 / t(bm)
df2 <- cbind(data$LT50d, df2)
names(df2)[names(df2) == 'data$LT50d'] <- 'LT50'
## Data frame 3 = CHC mass per body mass, CHC classes
df3 <- subset(data, select = c(89:93))/ t(bm)
df3 <- cbind(data$LT50d, df3)
names(df3)[names(df3) == 'data$LT50d'] <- 'LT50'
#

## BLOCK 1: PLSR models (byTube) -----------------------------------------------------------------------------------------

## PLSR data frames
dfp1 <- subset(data, select = c(LT50d,Avg,pL,pAke,pMo,pDi,pTri,wChain))
dfp2 <- subset(data, select = c(LT50d,Avg,SumL,SumAke,SumMo,SumDi,SumTri,wChain)) # raw mass estimates
dfp3 <- df3

#dfp3 <- subset(data, select = c(LT50d,SumL,SumAke,SumMo,SumDi,SumTri,wChain)) # CHC mass per body mass

## Partial Least Squares Regression (PLSR)
##https://cran.r-project.org/web/packages/pls/vignettes/pls-manual.pdf
## Model 1 - using CHC class proportions
plsr1 <- plsr(LT50d ~ ., data = dfp1, scale = TRUE, validation = "CV") # partial least squares regression
plsr1$loading.weights # Loading weight of each covariate on each component
## Model 2 - using CHC class mass estimates
plsr2 <- plsr(LT50d ~ ., data = dfp2, scale = TRUE, validation = "CV") # partial least squares regression
plsr2$loading.weights # Loading weight of each covariate on each component
## Model 3 - using CHC mass per body mass, summed for each class
plsr3 <- plsr(LT50 ~ ., data = dfp3, scale = TRUE, validation = "CV") # partial least squares regression
plsr3$loading.weights # Loading weight of each covariate on each component

## BLOCK 2: PLSR model RMSEP graph -----------------------------------------------------------------------------------------

## RMSEP of model 1 vs. model 2
ncomps <- c(1:7) # number of components
RM1 <- c(2.048,1.783,1.797,1.817,1.860,1.878,1.909) # Cross-validated RMSEP values for PLSR 1 (proportions)
RM2 <- c(2.113,1.945,1.894,1.871,1.841,1.895,1.779) # Cross-validated RMSEP values for PLSR 2 (mass estimates)
df <- data.frame(ncomps,RM1,RM2)
plot(ncomps,RM2,type = "b",xlab="Num. of components",ylab="RMSEP",pch=16,col="gray") # Plot with one model data points and line
lines(ncomps,RM1,type = "b",pch=16) # Add other model to same graph
legend("topright", legend = c("Model 1", "Model 2"), # Add legend explaining two models
       pch = c(16,16),col=c("black","gray"))

## BLOCK 3: PLSR loading weights to compare to random forest selection ---------------------------------------------------------

## Partial least squares regression
m1 <- plsr(LT50 ~ ., data = df1, scale = TRUE, validation = "CV") # partial least squares regression
## Data frame for forest plot showing all loading weights in a component
cov <- colnames(data)[17:88]
lws <- as.numeric(m1$loading.weights[,1]) # loading weights for a component
dfp <- data.frame(cov,lws)
colnames(dfp) <- c("COV","LWS")
norder <- dfp[order(dfp$LWS),][,1] # order of COV by ascending LWS
dfp$COV <- factor(dfp$COV, levels = c(norder))

## Hard threshold for loading weights (w* = median(w)/interquartile(w))
#https://analyticalsciencejournals-onlinelibrary-wiley-com.libproxy.berkeley.edu/doi/10.1002/cem.3226
ml <- median(dfp$LWS) # median of loading weights for Comp 1
qr <- 0.03372973-0.15120447 # 75% minus 25% from quantile(dfp$LWS)
w <- ml/qr # w* = median(w)/interquartile range(w)\

## Manual forest plot (conf. int. plot)
gHR <- ggplot(dfp, aes(x=COV, y=LWS)) +
  theme_bw() +
  #theme_black() +
  #theme(plot.margin=unit(c(1,2,1,2),"cm")) +
  geom_hline(yintercept=0, lty=2) +
  #geom_errorbar(aes(ymin=LOW, ymax=UPP), color="darkgrey", width=0.25) +
  geom_point(stat="identity", shape=1, size=3) +
  #scale_y_continuous(position = "right") + 
  coord_flip() +
  ylab("") +
  xlab("")
#

## BLOCK 4: Random forest with individual CHCs to predict LT50d -------------------------------------------------------
# Split the data set into training and test sets
split <- sample.split(df2, SplitRatio = 0.7)
train <- subset(df2, split == "TRUE") # 70% of data
test <- subset(df2, split == "FALSE") # 30% of data

set.seed(120)  # Setting seed. IIRC this is preferred before randomization in general
reg_RF = randomForest(x = train[-1], # remove the LT50d column
                      y = train$LT50, # because LT50d will be the Y we try to predict
                      ntree = 1000, type = "regression")
y_pred = predict(reg_RF, newdata = test[-1]) # also remove species column from test dataset
confusion_mtx = table(test[,1], y_pred)
#importance(reg_RF) # showing only top 30 variables, I think
varImpPlot(reg_RF)z
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
## BLOCK 5: Spearman correlating individual CHCs vs. LT50 -------------------------------------------------------------------
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
  y <- df1$LT50 # same Y variable each time
  x <- df1[,i+1] # current CHC column
  p <- as.numeric(cor.test(x,y,method="spearman",exact=FALSE)[3]) # p-value of current spearman test
  c <- as.numeric(cor.test(x,y,method="spearman",exact=FALSE)[4]) # p-value of current spearman test
  plist <- c(plist,p)
  clist <- c(clist,c)
}

# Spearman data frame using df1 CHC values (CHC mass)
sdf1 <- tibble(colnames(df1)[-1],p.adjust(plist,method="BH"),clist)
colnames(sdf1) <- c("CHC","p-adjusted","coef")
sdf1 <- sdf1[order(sdf1$`p-adjusted`),] # lowest value is most significant

# Spearman data frame using df2 CHC values (CHC mass per body mass)
sdf2 <- tibble(colnames(df2),p.adjust(plist,method="BH"),clist)
colnames(sdf2) <- c("CHC","p-adjusted","coef")
sdf2 <- sdf2[order(sdf2$`p-adjusted`),] # lowest value is most significant

## BLOCK 6: CHC stacked bar plot ----------------------------------------------------------------------------------------
## CHC class proportions df 
dfb <- read.csv("CHC_classes_StackedBar.csv")
dfb$NEST <- factor(dfb$NEST, levels = c("UK","DA","AB","LP","MT","LS","LH","SW"))
dfb$CLASS <- factor(dfb$CLASS, levels = c("Tri","Di","Mono","Alke","Alka"))
## Stacked bar plot - TOTAL M:L
sb1 <- ggplot(dfb, aes(fill=CLASS, y=SUM, x=NEST)) + 
  #theme_minimal() +
  theme(plot.title = element_text(size = 13, hjust = 0.5)) + # Change the title font
  labs(title="CHC class mass estimates") +
  geom_bar(position="dodge", stat="identity") +
  #geom_vline(aes(xintercept=3.5),linetype="dashed",size=0.5) +
  #scale_fill_manual(values=c("#CA2E55","#341566"),name="Class",labels=c("Total Methyl alkanes","Total n-Alkanes")) +
  #scale_fill_manual(values=c("#9E2A2B","#042A2B"),name="Class",labels=c("Total Methyl alkanes","Total n-alkanes")) +
  #scale_fill_manual(values=c("#341566","#51E18D","#6A83F1","#ECAA46","#9E2A2B")) +
  scale_fill_manual(values=c("#9E2A2B","#ECAA46","#6A83F1","#51E18D","darkorchid4")) +
  #theme(legend.position="none") +
  xlab("Nest") +
  ylab("Mass estimate (ng) on average profile")
#
