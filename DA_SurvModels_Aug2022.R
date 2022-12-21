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
library(GGally) # for ggpairs correlation plots
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

## BLOCK 2: Linear or logistic regression (Data = byTube) -------------------------------------------------------------------
m1 <- lm(LT50d ~ Avg, data = data)
m2 <- lm(LT50d ~ Avg + SumCHC, data = data)
m3 <- lm(LT50d ~ Avg + wChain, data = data)
m4 <- lm(LT50d ~ Avg + pL + pAke + pMo + pDi + pTri + wChain, data = data)

## BLOCK 3: PCR and PLSR models (byTube) -----------------------------------------------------------------------------------------

## PCR (Principal component regression) using pls package
##https://www.r-bloggers.com/2016/07/performing-principal-components-regression-pcr-in-r/
dfp1 <- subset(data, select = c(LT50d,Avg,pL,pAke,pMo,pDi,pTri,wChain))
dfp2 <- subset(data, select = c(LT50d,Avg,SumL,SumAke,SumMo,SumDi,SumTri,wChain))
dfp3 <- subset(data, select = c(LT50d,Avg,pL,pAke,pMo,pDi,pTri,Chain))
dfp4 <- subset(data, select = c(LT50d,Avg,SumL,SumAke,SumMo,SumDi,SumTri,Chain))
dfp5 <- subset(data, select = c(LT50d,Avg,pL,pAke,pMo,pDi,pTri,
                                L.chain,Ake.chain,Mo.chain,Di.chain,Tri.chain))
dfp6 <- subset(data, select = c(LT50d,L.chain,Ake.chain,Mo.chain,Di.chain,Tri.chain))
dfp7 <- subset(data, select = c(LT50d,pChain.L,pChain.Mo,pChain.Di,pChain.Tri))
dfp8 <- subset(data, select = c(LT50d,pL,pAke,pMo,pDi,pTri,wChain))
dfp9 <- subset(data, select = c(LT50d,19:90)) # all CHCs individually

## Partial Least Squares Regression (PLSR)
##https://cran.r-project.org/web/packages/pls/vignettes/pls-manual.pdf
plsr1 <- plsr(LT50d ~ ., data = dfp1, scale = TRUE, validation = "CV") # partial least squares regression
validationplot(plsr1, val.type = "R2") # according to R2 fit, only two PCA components are needed to make best model fit
plot(RMSEP(plsr1), legendpos = "topright") # How many components are needed to minimize RMSEP
plot(plsr1, ncomp = 1, asp = 1, line = TRUE) # Cross-validated predictions, using x num of components. asp = aspect ratio.
ncomp.permut <- selectNcomp(plsr1, method = "randomization", plot = TRUE)
plot(plsr1, plottype = "correlation")
biplot(plsr1, which = "loadings")
plsr1$loading.weights # Loading weight of each covariate on each component

## BLOCK 4: PLSR model comparison graphs -----------------------------------------------------------------------------------------

## RMSEP of model 1 vs. model 2
ncomps <- c(1:7) # number of components
RM1 <- c(2.048,1.783,1.797,1.817,1.860,1.878,1.909) # Cross-validated RMSEP values for PLSR 1 (proportions)
RM2 <- c(2.113,1.945,1.894,1.871,1.841,1.895,1.779) # Cross-validated RMSEP values for PLSR 2 (mass estimates)
df <- data.frame(ncomps,RM1,RM2)
plot(ncomps,RM2,type = "b",xlab="Num. of components",ylab="RMSEP",pch=16,col="gray") # Plot with one model data points and line
lines(ncomps,RM1,type = "b",pch=16) # Add other model to same graph
legend("topright", legend = c("Model 1", "Model 2"), # Add legend explaining two models
       pch = c(16,16),col=c("black","gray"))

## Data frame for Prediction plots
RV <- c(dfp1$LT50d,dfp1$LT50d) # Actual LT50
m1 <- data.frame(list(plsr1$validation$pred))[,2] # PLSR model 1 predictions, using 2 components
m2 <- data.frame(list(plsr2$validation$pred))[,2] # PLSR model 2 predictions, using 2 components
PV <- c(m2,m1)# Predicted values (from both models)
ID <- c(rep("Model 2",80),rep("Model 1",80))
#df <- data.frame(RV,PV,ID)
df <- data.frame(RV,m1,m2)

## Residuals for Predicted vs. Actual LT50 
#resid1 <- with(df, sum(m1-RV)^2)
#resid2 <- with(df, sum(m2-RV)^2)
#r1 <- resid(lm(RV~m1, data = df))

## Predicted vs. Actual LT50d
gs0 <- ggscatter(df, x = "PV", y = "RV",
                 add = "reg.line", conf.int = FALSE, 
                 cor.coef = FALSE, cor.method = "spearman",
                 xlab = "Predicted LT50", ylab = "Actual LT50",
                 color = "ID", palette = c("black","gray"),
                 xlim = c(0,16.5), ylim = c(0,16.5)) +
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))
  geom_abline(intercept = 0, slope = 1, linetype = 2) # 1:1 line, since this is Observed vs. Expected


## Predicted vs. Actual LT50d
gs1 <- ggscatter(df, x = "m1", y = "RV",
                 add = "reg.line", #conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "spearman",
                 xlab = "Predicted LT50", ylab = "Actual LT50",
                 xlim = c(0,16.5), ylim = c(0,16.5)) +
                 #color = "ID", palette = c(primary,cLH,cLS,cSW)) +
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))
  geom_abline(intercept = 0, slope = 1, linetype = 2) # 1:1 line, since this is Observed vs. Expected

## BLOCK 5: PLSR component lws vs. LT50d -----------------------------------------------------------------------------------------

## Data frame for loading weight correlations
df <- data.frame(list(plsr2$validation$pred))
Comp1 <- plsr2$scores[,1]
Comp2 <- plsr2$scores[,2]
RV <- dfp2$LT50d
#Res <- df$
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

## BLOCK 6: Tables (latex, formattable) -----------------------------------------------------------------------------------------
## LaTex tables
xdf <- read.csv("Table1_PLSR_Summary.csv")
print(xtable(xdf),
      include.rownames = FALSE,
      floating = TRUE, 
      latex.environments = "center",
      #floating.environment = "sidewaystable"
)
data(oliveoil)
mod <- plsr(sensory ~ chemical, data = oliveoil)
plot(mod, plottype = "biplot")

## Formattable example: https://cran.r-project.org/web/packages/formattable/vignettes/formattable-data-frame.html
#products <- data.frame(id = 1:5, 
#                       price = c(10, 15, 12, 8, 9),
#                       rating = c(5, 4, 4, 3, 4),
#                       market_share = percent(c(0.1, 0.12, 0.05, 0.03, 0.14)),
#                       revenue = accounting(c(55000, 36400, 12000, -25000, 98100)),
#                       profit = accounting(c(25300, 11500, -8200, -46000, 65000)))
#formattable(products, list(
#  price = color_tile("transparent", "lightpink"),
#  rating = color_bar("lightgreen"),
#  market_share = color_bar("lightblue")))
## Format table with LM and PLSR outputs
#xdf <- read.csv("ModelTable_lm_plsr.csv")
#psig <- formatter("span", style = x ~ style(color = ifelse(x < 0.05, "red", "black")))
#formattable(xdf, list(
#  p.value = psig,
#  COMP.1 = color_tile("transparent","lightpink"),
#  COMP.2 = color_tile("transparent","lightpink")
#))
