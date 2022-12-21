***

## I: New solutions to the multicollinearity problem 
***

The JEB reveiwers had some good points, but the strongest criticisms are aimed at using alkane ratios in our cox regression model, as well as using individual-level survival data (i.e. Time-to-death of each individual ant) when our covariates are at the level of the replicate group or nest. These were concerns of mine, too, as we built the paper around these cox regression models. Ultimately, this structure of our data, and our methods of analysis, were attempts to solve the problem of multicollinearity in our CHC data. CHCs are not independent variables, but if we want to know which CHCs or classes explained survival, we needed to treat them as separate variables in a regression model, which lead to impossibly large coefficients that cannot be trusted/interpreted. Sometimes, when covariates are so correlated with each other, regression models can't even output coefficients, and just give errors or "NA" instead. The ratio was an attempt to circumvent this by condensing all of the correlated CHC data into one or two covariates, but it is indeed an improper way to test for opposing effects from CHC class. Setting them up as a ratio means the compound classes would always have opposing effects (numerator vs. denominator) if they had any effects at all. 

So, we need a new strategy for dealing with multicollinearity, if we do want to examine CHC proportions as covariates that explain survival during desiccation. I read some relevant chapters from ["Regression Model Strategies"](https://link.springer.com/book/10.1007/978-3-319-19425-7), and came across the solution many of us have heard of before: putting the data in a PCA model, and using the principal components as the model covariates instead. Previously, I had avoided this because the interpretation of a PCA and its components is limited. A PCA allows you to examine the dimensionality of the X values, and if you show that colonies are different in certain X values, then you could suggest that maybe this difference also explained the survival differences (Y values). However, there is a version of principal component regression that not only transforms X data into different dimensions, but also estimates how those dimensions relate to the response variable (Y values).

This method is called **"Partial Least Squares Regression" (PLSR)**. It is called this because it partially involves using least squares in its formula to determine a line of best fit. An odd name, since I don't think it explains what makes it special. Here is a good description from [this website](http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/152-principal-component-and-partial-least-squares-regression-essentials/): 

> *"A possible drawback of PCR [principal component regression] is that we have no guarantee that the selected principal components are associated with the outcome. Here, the selection of the principal components to incorporate in the model is not supervised by the outcome variable. An alternative to PCR is the Partial Least Squares (PLS) regression, which identifies new principal components that not only summarizes the original predictors, but also that are related to the outcome. These components are then used to fit the regression model. So, compared to PCR, PLS uses a dimension reduction strategy that is supervised by the outcome."*

Therefore, a PLS regression (PLSR) is a better tool if we want to claim that certain components or covariates drove survival in our desiccation assays. If you want to fully understand this method and its use in ecological research like our own, I highly suggest reading [this paper](https://onlinelibrary-wiley-com.libproxy.berkeley.edu/doi/10.1111/j.1600-0706.2008.16881.x) (cited below), as it explains why PLS is useful and how to interpret its results. For now, I will jump into performing the analysis in our data, and how to interpret the results. 

*Carrascal, L. M., Galván, I., & Gordo, O. (2009). Partial least squares regression as an alternative to current regression methods used in ecology. Oikos, 118(5), 681-690.*

***
## II: Data frame and PLSR summary 
***

```{r eval=TRUE}
## Here is the structure or scale of our data, modified for use in PLS instead of Cox Regression:
dfp <- subset(data, select = c(LT50d,Avg,pL,pMo,pDi,pTri,wChain))
head(dfp)
```
Note that Time-to-death of each ant is no longer the survival variable. Instead, we are using the LT50 of each test tube in the drierite treatment, because I think we shouldn't use individual-level survival data if none of our covariates are also at the individual-level. Also, remember we only use the drierite data in these models that include the body size, because the drierite tubes were the only ones we saved and measured body size. Thus, LT50 and Avg (body size) are unique to each row. The remaining variables are CHC metrics, which are the same for each nest ID. A list of abbreviations:

- pl = % of n-alkanes on profile
- pMo = % of mono-methyl alkanes on profile
- pDi = % of di-methyl alkanes on profile
- pTri = % of tri-methyl alkanes on profile
- wChain = Weighted average chain length of profile

All of the CHC and body size data is the same as what we were using for the Cox regression, but no longer as ratios, and no longer divided into selected vs. total M:L. The proportions of each class for each nest ID come from Jan's supplemental data of the 2018 publication, which contains the nanogram estimates of each identified CHC from our ant nests. We no longer only examine a selected subset of the CHCs, but we can, if we want to. Inputting this data into a PLS looks like this:

```{r eval=TRUE}
## Partial least squares regression (formula + summary)
## https://cran.r-project.org/web/packages/pls/vignettes/pls-manual.pdf
plsr_model <- plsr(LT50d ~ ., data = dfp, scale = TRUE, validation = "CV")
## LT50 vs. all covariates in dfp. Scale = TRUE scales the data, similar to the z-scoring we did before. Validation = CV is for cross-validation. Other validation options seem more specialized. 
summary(plsr_model)
```
Just like a principal components analysis, the covariates are replaced with components that are transformations of the covariates you provided. Ideally, the majority of the variance in your data can be explained with a few of these components. So the first step to interpretation is to decide how many components are needed to explain our survival data. Remember, PLSR allows us to not only see how much of the X variance is explained (body size, CHC proportions, etc.) but also how much of the Y variance is explained (LT50 from all drierite tubes) by these components. The variance explained can be seen in the summary above, but I will summarize it in the table below:

<center>
![](H:/My Drive/R_folder/1 - NeilNSF/1 - Publication/Resubmitting/Table1_PLSRSummary.jpg)

</center>

PC1 tells us the most about LT50, and the second component also adds some detail, but beyond that the help becomes negligible. PC2 is useful if we want to understand the dimensionality of the X variables (which is what usual PCA biplots show). RMSEP stands for "Root Mean Square Error of Prediction", and the decision of how many components to select should be judged by how minimized we can make this error prediction. Because 1.81 is as low as we can go, and we can achieve this with just the first two components, it is safe to just use these two components for our analysis and interpretation.

Each of these components is made from transformations of the covariates, and some of the covariates had more "weight" in these transformations. This part is hard to explain (some background provided in [this video](https://youtu.be/Vf7doatc2rA)), but to put it simply, if PC1 explains most of the Y variance, we can look at the "loading weights" from each covariate for PC1 in order to see which ones correlate most with PC1. Perhaps this will be easier to think of if we look at a regular PCA biplot using our data:

<center>
![](H:/My Drive/R_folder/1 - NeilNSF/1 - Publication/Resubmitting/XBiplot.jpeg)

</center>

The black numbers are the X scores (one for each sample and each component) and the red arrows are the loadings, or direction vectors. The direction vectors suggest that placement on the x axis (Comp 1) is largely influenced by body size and % of di-methyl alkanes, while the other covariates determine the y axis (Comp 2) placement of the scores. Without further information, we can use this to suggest that body size and % of di-methyl alkanes oppose each other, so if body size helps survival, then di-methyl alkanes might do the opposite. However, like the quote I mentioned before, there is no guarantee that the dimensionality of these components are associated with the outcome variable (survival). 

Because this is a PLSR, though, we know Comp 1 explains LT50, and we can see the loading weights of each covariate in Comp 1 to claim which covariates most likely positively or negatively influenced survival. Making claims like this is similar to how we interpret the covariates or predictors of a multiple linear regression, so I'll show both the PLSR loading weights and the summary of a multiple regression below (using the basic `lm( )` function in R).

***
## III: PLSR results and interpretation 
***

<center>
![](H:/My Drive/R_folder/1 - NeilNSF/1 - Publication/Resubmitting/Table2_LM_PLSR.jpg)

</center>

The results of these two analyses are similar, but the multiple regression can't be reasonably interpreted. More specifically, its coefficients are ridiculously large, and the % of tri-methyl alkanes can't even be considered because it is so correlated with the other covariates. The fact that we still have significant values is actually a good sign, because usually inflated coefficients also leads to inflated standard errors and overlapping confidence intervals. The PLSR, being a version of PCA, removes all problems of multicollineairity because the covariates are transformed into an orthogonal matrix, and the loading weights of each covariate in each component does not require a p-value to be considered relevant to interpretation. 

For Comp 1, body size has the strongest positive effect and % of di-methyl alkanes has the strongest negative effect. We saw this on the biplot. Since Comp 1 explains the Y variance the best, higher body size is associated with increased survival, while a higher % of di-methyl alkanes is associated with decreased survival. This makes sense when we look at Comp 1 vs. LT50:

<center>
![](H:/My Drive/R_folder/1 - NeilNSF/1 - Publication/Resubmitting/PC1_vs_LT50d.jpg)

</center>

So Comp 1 really aligns with survival during desiccation, and we know which covariates define Comp 1. The remaining covariates help explain a different dimension of the data. Instead of correlating well with LT50, Comp 2 seems to separate the northern and southern nests in terms of their remaining CHC properties (% n-alkanes, mono-methyl, tri-methyl, and chain lengths):

<center>
![](H:/My Drive/R_folder/1 - NeilNSF/1 - Publication/Resubmitting/PC2_vs_LT50d.jpg)

</center>

This is as far as my interpretation and understanding go, right now. Perhaps this new analysis invites new investigations into the CHC profiles. For instance, does it make sense to consider di-methyl alkanes has being worse at water proofing than the other methyl alkanes? Are the di-methyl alkanes in that weak-zone of not only having methyl branches, but also having small chain lengths? Last I remember, they were around the same size as the tri-methyl alkanes on the profile. 

We can also consider other covariates to model. I chose the categories above, as it seemed useful to consider more than just methylated vs. not-methylated alkanes, like we were stuck doing with the ratios before. By effectively removing the problem of multicollinearity, we can freely investigate which features of the CHC profile might have more of an influence than other features. If we are less restricted in what covariates we want to explore, we might need to consider methods for choosing models (like AIC scoring) instead of sticking to one model based on one hypothesis of interest. 

***
