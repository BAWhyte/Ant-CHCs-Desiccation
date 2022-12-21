########################################################
### Partial least squares regression for survival analysis ###
# created by Brian Whyte, Sep 2021
# 
###  - - - - - - - - - - - -   ###
## BLOCK 0: Set up + clean up -------------------------------------------------------

## Housekeeping
cat("\014") # Clear console
rm(list=ls()) # Remove all variables

## Libraries
library(ggplot2) # plotting
library(ggpubr) # for ggscatter plots
library(survival) # for survival functions
library(dplyr) # data transforming
library(pls) # for principle component regression
library(stargazer) # for automatic latex outputs of lm objects
library(xtable) # for latex table outputs
library(ggbiplot) # for ggplot type PCA biplots

## Set WD + load data
setwd("Your/Path/")
data <- read.csv("ExampleData.csv")

## Colony colors
primary <- "#FFDF66" # Main (large) supercolony
cLH <- "#0582CA" # Lake Hodges
cLS <- "mediumorchid1" # Lake Skinner
cSW <- "lightgreen" # SweetWater

## BLOCK 1: Linear and PLSR models (byTube) -----------------------------------------------------------------------------------------
## Linear regression (not trustworthy, because of collinearity)
lm1 <- lm(LT50d ~ Avg + pL + pAke + pMo + pDi + pTri + wChain, data = data)

## Partial Least Squares Regression (PLSR)
##https://cran.r-project.org/web/packages/pls/vignettes/pls-manual.pdf
dfp1 <- subset(data, select = c(LT50d,Avg,pL,pAke,pMo,pDi,pTri,wChain)) # isolate data frame of only selected covariates (+ response variable)
plsr1 <- plsr(LT50d ~ ., data = dfp1, scale = TRUE, validation = "CV") # partial least squares regression
validationplot(plsr1, val.type = "R2") # according to R2 fit, only two PCA components are needed to make best model fit
plot(RMSEP(plsr1), legendpos = "topright") # How many components are needed to minimize RMSEP
plot(plsr1, ncomp = 1, asp = 1, line = TRUE) # Cross-validated predictions, using x num of components. asp = aspect ratio.
ncomp.permut <- selectNcomp(plsr1, method = "randomization", plot = TRUE)
plot(plsr1, plottype = "correlation")
biplot(plsr1, which = "loadings")
plsr1$loading.weights # Loading weight of each covariate on each component

## BLOCK 2: PLSR component lws vs. LT50d -----------------------------------------------------------------------------------------

## Data frame for loading weight correlations
df <- data.frame(list(plsr2$validation$pred))
Comp1 <- plsr2$scores[,1]
Comp2 <- plsr2$scores[,2]
RV <- dfp2$LT50d
ID <- data$Super
Name <- data$NestID
df <- data.frame(cbind(ID,Name,RV,Comp1,Comp2,df))

## Spearman correlation of PLSR components vs. LT50d
gs1 <- ggscatter(df, x = "Comp1", y = "RV",
                 label = "Name", repel = TRUE,
                 #add = "reg.line", #conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "spearman",
                 color = "ID", palette = c(primary,cLH,cLS,cSW),
                 xlab = "PSLR Component 1", ylab = "Mean LT50 per replicate tube (Drierite only)") +
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_smooth(method=lm, se=FALSE, color = "grey70", linetype = 2)
  #theme_black()
gs2 <- ggscatter(df, x = "Comp2", y = "RV", 
                 label = "Name", repel = TRUE,
                 #add = "reg.line", #conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "spearman",
                 color = "ID", palette = c(primary,cLH,cLS,cSW),
                 xlab = "PSLR Component 2", ylab = "Mean LT50 per replicate tube (Drierite only)") +
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_smooth(method=lm, se=FALSE, color = "grey70", linetype = 2) +
  theme_black()

## BLOCK 3: Tables (latex, formattable) -----------------------------------------------------------------------------------------
## LaTex tables
xdf <- read.csv("PLSR_TableData.csv")
print(xtable(xdf),
      include.rownames = FALSE,
      floating = TRUE, 
      latex.environments = "center",
      #floating.environment = "sidewaystable"
)
