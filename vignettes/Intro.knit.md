---
title: "Introduction to fitting statistical models to collections of enhancers"
author: "Noah Dukler"
date: "2018-04-02"
output: pdf_document
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

\newcommand{\vect}[1]{\boldsymbol{\mathbf{#1}}}

## Statistical models for super-enhancers

Activity model is $A(\vect{x}) = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \cdots + \beta_n
x_n$.

Linear model is $R(\vect{x}) = e^{A(\vect{x})} + \epsilon$.

Exponential model is $R(\vect x) = e^{A(\vect x)} + \epsilon$.

Logistic model is $R(\vect x)=1/\left(1+e^{-A(\vect{x})}\right)+\epsilon$.

## Fitting your first model

To start fitting statistical models to your enhancer data you need to format your data in a particular way. First you need to place information about which enhancers are present in each experiment in a data.frame where each row is an experiment and each column is an enhancer element. Here's a simple example using the data from Shin et al. (2015):


```r
library(superEnhancerModelR)
data("wap")
design = wap[,c(3:5)]
head(design)
```

```
##   E1 E2 E3
## 1  1  1  1
## 2  1  1  1
## 3  1  1  1
## 4  1  1  1
## 5  1  1  1
## 6  1  1  1
```

Then we need a vector of normalized expression values, with each element of the vector corresponding to one row of the previous data frame.


```r
## Create "expression data""
expression = wap[,2]
print(head(expression))
```

```
## [1] 141.60 110.68  47.72 114.84 130.00  62.14
```

Once we have both expression values and a design matrix, we can construct an enhancerDataObject. This object holds the raw data, the design of each experiment, and detail of the statistical model that you want to fit. You can specify the link function: additive, exponential, or logistic (from Dukler et al., 2017). You can also specify an error model: Gaussian, or log-normal. Upper and lower bounds for all parameters can also be specified although we will leave them at the defaults for now. We specify a model with a logistic link function, log-normal error, and non-interacting enhancer loci as follows:


```r
## Create activity function
actFun = formula(~E1+E2+E3)
independent.loci = enhancerDataObject(expression,design,actFun,errorModel="lognormal",linkFunction="logistic")
```

This model can then be fit using a differential evolutionary algorithm and refined with gradient descent (note that this step could be run with more iterations.):

```r
independent.loci=optimDE(independent.loci,maxit=500,refine=TRUE)
```

```
## Iteration: 100 bestvalit: 180.829285 bestmemit:   -5.258533  109.135808  454.528025  757.333183   26.275068    1.153662
## Iteration: 200 bestvalit: 178.071923 bestmemit:   -4.954458   85.512336  -80.179712  596.583499   25.970310    1.053851
## Iteration: 300 bestvalit: 175.023113 bestmemit:   -5.579245  124.633513 -120.181650  556.089250   33.027127    0.973075
## Iteration: 400 bestvalit: 169.077595 bestmemit:   -5.607358    4.267754    0.911145  924.327670   38.987021    0.905565
## Iteration: 500 bestvalit: 165.517263 bestmemit:   -5.826014    4.144612    0.163777  351.359569   41.652620    0.787104
```

Once you have a fitted model you can plot the data and the predicted model


```r
plotModel(independent.loci)
```

![](C:\Users\Noah Dukler\AppData\Local\Temp\RtmpEvh8G3\preview-3e10763e1d26.dir\Intro_files/figure-latex/unnamed-chunk-5-1.pdf)<!-- --> 

You can also view the model residuals

```r
plotResiduals(independent.loci)
```

![](C:\Users\Noah Dukler\AppData\Local\Temp\RtmpEvh8G3\preview-3e10763e1d26.dir\Intro_files/figure-latex/unnamed-chunk-6-1.pdf)<!-- --> 


Say we want to test whether E1 and E2 interact. Then we augment the independent model with an interaction term (E1*E2), and fit it. We can view the fitted interaction model in the same way we viewed the independent loci model.


```r
actFun2 = formula(~E1+E2+E1*E2)
dependent.loci = enhancerDataObject(expression,design,actFun2,errorModel="lognormal",linkFunction="logistic")
dependent.loci=optimDE(dependent.loci,maxit=300,refine=TRUE)
```

```
## Iteration: 100 bestvalit: 191.791710 bestmemit:   -2.742207  266.117690  559.678306  228.809560   22.520833    1.700698
## Iteration: 200 bestvalit: 190.054429 bestmemit:   -3.182886    4.600900  320.117350 -124.873427   30.323660    1.396322
## Iteration: 300 bestvalit: 190.029072 bestmemit:   -3.208540    4.164233    9.298622  151.642480   31.386067    1.405659
```

```r
plotModel(dependent.loci)
```

![](C:\Users\Noah Dukler\AppData\Local\Temp\RtmpEvh8G3\preview-3e10763e1d26.dir\Intro_files/figure-latex/unnamed-chunk-7-1.pdf)<!-- --> 
Then we can compare their respective Bayesian information criterion (BIC). Lower BIC indicates a better model.

```r
paste("The BIC for the interacting loci model is:",bic(dependent.loci))
```

```
## [1] "The BIC for the interacting loci model is: 402.339569509411"
```

```r
paste("The BIC for the independent loci model is:",bic(independent.loci))
```

```
## [1] "The BIC for the independent loci model is: 352.075113405854"
```

The significance of a difference in Bayes factor can be interpreted using this table.

$\Delta BIC$   Evidence against higher BIC        
-------------  -----------------------------------
0 to 2         Not worth more than a bare mention 
2 to 6         Positive                           
6 to 10        Strong                             
>10            Very Strong                        

If you want p-values for individual coefficients, run the function computeParamPvals and it will calculate the significance of the each parameter via the likelihood-ratio test.  


```r
param.stats=computeParamPvals(independent.loci,maxit=500)
```

```
## Fitting reduced model without term (Intercept)
## Iteration: 100 bestvalit: 195.458917 bestmemit:   -1.579047  528.962721  944.731343  441.407235   21.272770    1.552989
## Iteration: 200 bestvalit: 193.714687 bestmemit:   -1.824591    1.128457  515.635597  223.521097   26.016840    1.538323
## Iteration: 300 bestvalit: 193.686872 bestmemit:   -1.689803    0.822383  148.065658 -388.206550   26.316917    1.536791
## Iteration: 400 bestvalit: 193.686871 bestmemit:   -1.689333    0.821457  334.131243 -571.267218   26.314017    1.536917
## Iteration: 500 bestvalit: 193.686871 bestmemit:   -1.689337    0.821463   88.997187  323.131996   26.314027    1.536917
## Fitting reduced model without term E1
## Iteration: 100 bestvalit: 183.653270 bestmemit:   -3.618727  505.830586  120.032639  835.947986   29.854076    1.145788
## Iteration: 200 bestvalit: 176.811165 bestmemit:   -4.086567    2.425482  465.982516 -585.863544   43.580552    1.005264
## Iteration: 300 bestvalit: 175.277874 bestmemit:   -4.454467    2.469363    4.383800  385.518311   62.766680    0.982511
## Iteration: 400 bestvalit: 175.272238 bestmemit:   -4.479723    2.523446    4.383874  295.517234   62.803255    0.981221
## Iteration: 500 bestvalit: 175.272217 bestmemit:   -4.478989    2.523830    4.382166  135.573161   62.758184    0.980828
## Fitting reduced model without term E2
## Iteration: 100 bestvalit: 180.887923 bestmemit:   -5.066502  138.258017  305.416340  772.438951   25.673095    1.124216
## Iteration: 200 bestvalit: 166.036720 bestmemit:   -5.541676    3.979963  346.234214  646.267289   40.946697    0.809900
## Iteration: 300 bestvalit: 160.293262 bestmemit:   -6.417145    4.171473    5.963560 -312.045596   64.373913    0.690301
## Iteration: 400 bestvalit: 160.151585 bestmemit:   -6.254631    4.087064    5.753128 -102.724805   62.661313    0.678782
## Iteration: 500 bestvalit: 160.151400 bestmemit:   -6.250237    4.079026    5.748043 -794.688350   62.723630    0.678284
## Fitting reduced model without term E3
## Iteration: 100 bestvalit: 190.352499 bestmemit:   -3.278543  785.952193   39.478263   82.364824   26.475639    1.413376
## Iteration: 200 bestvalit: 190.081733 bestmemit:   -3.071253    4.422483   69.201241 -868.422221   28.549418    1.406727
## Iteration: 300 bestvalit: 190.029084 bestmemit:   -3.207766    4.162312    9.549384  315.224837   31.409873    1.405436
## Iteration: 400 bestvalit: 190.029061 bestmemit:   -3.206909    4.162939    8.634904  756.946148   31.387396    1.405741
## Iteration: 500 bestvalit: 190.029061 bestmemit:   -3.206983    4.162948    8.639737   87.456775   31.388901    1.405740
```

```r
param.stats$pVals
```

```
##     parameter        ll        D        p.val
## 1 (Intercept) -193.6869 57.58006 1.649767e-14
## 2          E1 -175.2722 20.75075 2.731645e-06
## 3          E2 -160.1514 -9.49088 0.000000e+00
## 4          E3 -190.0291 50.26444 6.846913e-13
```

