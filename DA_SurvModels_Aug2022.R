### ANT CHCS AND DESICCATION ###
# created by Brian Whyte, Sep 2020
# see repository for up-to-date README
###  - - - - - - - - - - - -   ###
## BLOCK 0: Set up + clean up -------------------------------------------------------

## Housekeeping
cat("\014") # Clear console
rm(list=ls()) # Remove all variables

## Libraries
library(ggplot2) # plotting
library(GGally) # for ggpairs correlation plots
library(ggpubr) # for ggscatter plots
library(survival) # for survival functions
library(dplyr) # data transforming
library(pls) # for principle component regression
library(stargazer) # for automatic latex outputs of lm objects
library(xtable) # for latex table outputs
library(ggbiplot) # for ggplot type PCA biplots
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

## Set WD + load data
#setwd("F:/Users/Brian/Google Drive (ba.whyte@berkeley.edu)/R_folder/1 - NeilNSF/1 - Publication")
setwd("H:/My Drive/0 - R/1 - NeilNSF/1 - Publication/Resubmitting")
#setwd("/Volumes/GoogleDrive/My Drive/0 - R/1 - NeilNSF/1 - Publication/Resubmitting")
data <- read.csv("ModelData_byTube_Aug2022.csv")
#data <- read.csv("ModelData_PrelimPCA_byClass.csv")
#data1 <- read.csv("ModelData_byTube_July2022.csv")
#data2 <- read.csv("ModelData_byNest_Jun2022.csv")
#data2 <- subset(data2, select = c(NestID, Chain))
#df2 <- left_join(data1, data2, by = "NestID")
#write.csv(df2, "ModelData_byTube_Aug2022.csv")

## Colony colors
primary <- "#FFDF66" # Main (large) supercolony
cLH <- "#0582CA" # Lake Hodges
cLS <- "mediumorchid1" # Lake Skinner
cSW <- "lightgreen" # SweetWater

## BLOCK 1: Collinearity testing --------------------------------------------------------------------------------------------
# Collinearity plot (ggpairs) of chemical data
df1 <- subset(data, select = c(LT50d,Avg,pL,pMo,pDi,pTri,wChain))
df2 <- subset(data, select = c(LT50d,pL,pMo,pDi,pTri))
df3 <- subset(data, select = c(LT50d,Avg,SumCHC))
df4 <- subset(data, select = c(Avg,SumL,SumAke,SumMo,SumDi,SumTri,SumCHC))
#colnames(df1) <- c("Body.Size","% n-alka","% mono-me","% di-me","% tri-me","W.Chain")
#windows();
ggpairs(df4, lower = list(continuous = wrap("smooth")), upper = list(continous = wrap("cor", method = "spearman")))

# Spearman correlation X vs. Y
gs0 <- ggscatter(df1, x = "wChain", y = "LT50d", 
                 add = "reg.line", #conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "spearman")

## BLOCK 2: PLSR models (byTube) -----------------------------------------------------------------------------------------

## PCR (Principal component regression) using pls package
##https://www.r-bloggers.com/2016/07/performing-principal-components-regression-pcr-in-r/
dfp1 <- subset(data, select = c(LT50d,Avg,pL,pAke,pMo,pDi,pTri,wChain))
dfp2 <- subset(data, select = c(LT50d,Avg,SumL,SumAke,SumMo,SumDi,SumTri,wChain))
#dfp3 <- subset(data, select = c(LT50d,19:90)) # all CHCs individually
#dfp3 <- subset(data, select = c(LT50d,Avg,pL,pMethyl,wChain))
#dfp4 <- subset(data, select = c(LT50d,Avg,Reb_Alkanes,Reb_Malkanes,wChain))
#dfp <- data[c(17:88)] # All CHCs for preliminary PCA
#dfp <- data[c(3:10)] # All CHCs for preliminary PCA

## Preliminary PCA of all individual CHCs
#dfp.class <- data[,2]
#dfp.pca <- prcomp(dfp, center = TRUE, scale. = TRUE)
#gb <- ggbiplot(dfp.pca, obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = TRUE)

## Partial Least Squares Regression (PLSR)
##https://cran.r-project.org/web/packages/pls/vignettes/pls-manual.pdf
## Model 1 - using CHC class proportions
plsr1 <- plsr(LT50d ~ ., data = dfp1, scale = TRUE, validation = "CV") # partial least squares regression
validationplot(plsr1, val.type = "R2") # according to R2 fit, only two PCA components are needed to make best model fit
plot(RMSEP(plsr1), legendpos = "topright") # How many components are needed to minimize RMSEP
plot(plsr1, ncomp = 1, asp = 1, line = TRUE) # Cross-validated predictions, using x num of components. asp = aspect ratio.
ncomp.permut <- selectNcomp(plsr1, method = "randomization", plot = TRUE)
plot(plsr1, plottype = "correlation")
biplot(plsr1, which = "loadings")
plsr1$loading.weights # Loading weight of each covariate on each component
## Model 2 - using CHC class mass estimates
plsr2 <- plsr(LT50d ~ ., data = dfp1, scale = TRUE, validation = "CV") # partial least squares regression
validationplot(plsr2, val.type = "R2") # according to R2 fit, only two PCA components are needed to make best model fit
plot(RMSEP(plsr2), legendpos = "topright") # How many components are needed to minimize RMSEP
plot(plsr2, ncomp = 1, asp = 1, line = TRUE) # Cross-validated predictions, using x num of components. asp = aspect ratio.
ncomp.permut <- selectNcomp(plsr2, method = "randomization", plot = TRUE)
plot(plsr2, plottype = "correlation")
biplot(plsr2, which = "loadings")
plsr2$loading.weights # Loading weight of each covariate on each component

## BLOCK 3: PLSR model comparison graphs -----------------------------------------------------------------------------------------

## RMSEP of model 1 vs. model 2
ncomps <- c(1:7) # number of components
RM1 <- c(2.048,1.783,1.797,1.817,1.860,1.878,1.909) # Cross-validated RMSEP values for PLSR 1 (proportions)
RM2 <- c(2.113,1.945,1.894,1.871,1.841,1.895,1.779) # Cross-validated RMSEP values for PLSR 2 (mass estimates)
df <- data.frame(ncomps,RM1,RM2)
plot(ncomps,RM2,type = "b",xlab="Num. of components",ylab="RMSEP",pch=16,col="gray") # Plot with one model data points and line
lines(ncomps,RM1,type = "b",pch=16) # Add other model to same graph
legend("topright", legend = c("Model 1", "Model 2"), # Add legend explaining two models
       pch = c(16,16),col=c("black","gray"))

## BLOCK 4: PLSR component lws vs. LT50d -----------------------------------------------------------------------------------------

## Data frame for loading weight correlations
df <- data.frame(list(plsr2$validation$pred))
Comp1 <- plsr2$scores[,1]
Comp2 <- plsr2$scores[,2]
RV <- dfp2$LT50d
ID <- data$Super
Name <- data$NestID
df <- data.frame(cbind(ID,Name,RV,Comp1,Comp2,df))

## PLSR Comp 1 vs. LT50d
gs1 <- ggscatter(df, x = "Comp1", y = "RV",
                 label = "Name", repel = TRUE,
                 #add = "reg.line", #conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "spearman",
                 color = "ID", palette = c(primary,cLH,cLS,cSW),
                 xlab = "PSLR Component 1", ylab = "Mean LT50 per replicate tube (Drierite only)") +
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_smooth(method=lm, se=FALSE, color = "grey70", linetype = 2)
## PLSR Comp 2 vs. LT50d
gs2 <- ggscatter(df, x = "Comp2", y = "RV", 
                 label = "Name", repel = TRUE,
                 #add = "reg.line", #conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "spearman",
                 color = "ID", palette = c(primary,cLH,cLS,cSW),
                 xlab = "PSLR Component 2", ylab = "Mean LT50 per replicate tube (Drierite only)") +
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_smooth(method=lm, se=FALSE, color = "grey70", linetype = 2)

## BLOCK 5: CHC stacked bar plot ----------------------------------------------------------------------------------------
## CHC class proportions df 
dfb <- read.csv("CHC_classes_StackedBar.csv")
dfb$NEST <- factor(dfb$NEST, levels = c("UK","DA","AB","LP","MT","LS","LH","SW"))
dfb$CLASS <- factor(dfb$CLASS, levels = c("Tri","Di","Mono","Alke","Alka"))
## Stacked bar plot - TOTAL M:L
sb1 <- ggplot(dfb, aes(fill=CLASS, y=PERC, x=NEST)) + 
  #theme_minimal() +
  theme(plot.title = element_text(size = 13, hjust = 0.5)) + # Change the title font
  labs(title="CHC class proportions") +
  geom_bar(position="dodge", stat="identity") +
  #geom_vline(aes(xintercept=3.5),linetype="dashed",size=0.5) +
  #scale_fill_manual(values=c("#CA2E55","#341566"),name="Class",labels=c("Total Methyl alkanes","Total n-Alkanes")) +
  #scale_fill_manual(values=c("#9E2A2B","#042A2B"),name="Class",labels=c("Total Methyl alkanes","Total n-alkanes")) +
  #scale_fill_manual(values=c("#341566","#51E18D","#6A83F1","#ECAA46","#9E2A2B")) +
  scale_fill_manual(values=c("#9E2A2B","#ECAA46","#6A83F1","#51E18D","darkorchid4")) +
  #theme(legend.position="none") +
  xlab("Nest") +
  ylab("Relative proportions")
