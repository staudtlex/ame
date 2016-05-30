# ame 

**ame** computes simulation based (Krinsky-Robb) _average marginal effects_ and _first differences_ for various estimators in R. 

Main goals: compute marginal effects and standard errors for:
+ continuous variables (dy/dx)
+ dichotomous variables (first differences)
+ interactions of continuous variables (dy/dx)
+ interactions of categorial variables (first differences)
+ spline-functions (for generalized additive models)

## Installation
[![Build Status](https://travis-ci.org/staudtlex/ame.svg?branch=master)](https://travis-ci.org/staudtlex/ame)

To install `ame`, simply run

```R
devtools::install_github("staudtlex/ame")
```

## Examples
#### 1. Logit model, graduate school admission data
```R
library(ame)

# graduate school data
gsa <- read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
gsa$rank <- as.factor(gsa$rank)

# run glm logit model
glm_fit <- glm(admit ~ gre + gpa + rank, data = gsa, family = binomial(link = "logit"))

# compute average marginal effects
glm_ame <- ame(glm_fit, cont_vars = c("gre", "gpa"), fac_vars = c("rank"), nsim = 1000)
```

Use the `summary()` command to display the average marginal effects and first differences of the model variables:
```
Continuous variables: gre, gpa 
Factor variables:     rank 

             dydx  Std. Error z value  Pr(>|z|)    
gre    0.00043102  0.00020583  2.0941 0.0355246 *  
gpa    0.15578509  0.06248429  2.4932 0.0124070 *  
rank2 -0.15452167  0.07524027 -2.0537 0.0392037 *  
rank3 -0.28434508  0.07528385 -3.7770 0.0001556 ***
rank4 -0.31680193  0.08239552 -3.8449 0.0001182 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Note: For a factor variable f, dydx corresponds to the first difference E(Y|f_i) - E(Y|f_0)
```
#### 2. Linear model, mtcars-data
```R
# linear model using glm
lm1 <- glm(mpg ~ cyl * hp + wt, data = mtcars, family = gaussian)
summary(lm1)
```
```
Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept) 52.017520   4.916935  10.579 4.18e-11 ***
cyl         -2.742125   0.800228  -3.427  0.00197 ** 
hp          -0.163594   0.052122  -3.139  0.00408 ** 
wt          -3.119815   0.661322  -4.718 6.51e-05 ***
cyl:hp       0.018954   0.006645   2.852  0.00823 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
Compute average marginal effects and first differences.
```R
lm1_ame <- ame(lm2, cont_vars = c("cyl", "hp", "wt"), nsim = 1000, seed = 99)
summary(lm1_ame)
```
Display results:
```
Continuous variables: cyl, hp, wt 
Factor variables:      

         dydx Std. Error z value  Pr(>|z|)    
cyl  0.031662   0.602009  0.0526  0.938894    
hp  -0.046123   0.014881 -3.0993  0.001901 ** 
wt  -3.133662   0.650566 -4.8168 1.429e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Note: For a factor variable f, dydx corresponds to the first difference E(Y|f_i) - E(Y|f_0)
```

# Note
This version of `ame` supports `glm` objects with 
+ `gaussian` model family with `"identity"` link
+ `binomial` model family with `"logit"`/`"probit"` link
