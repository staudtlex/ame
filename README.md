# ame 

**ame** computes simulation based (Krinsky-Robb) _average marginal effects_ and _first differences_ for various estimators in R. 

Main goals: compute marginal effects and standard errors for:
+ continuous variables (dy/dx)
+ dichotomous variables (first differences)
+ interactions of continuous variables (dy/dx)
+ interactions of categorial variables (first differences)
+ spline-functions (for generalized additive models)

## Installation
To install `ame`, simply run

```R
devtools::install_github("staudtlex/ame")
```

## Examples
#### 1. Logit model with fake data
```R
library(ame)

# generate fake data
set.seed(100)
y <- rbinom(100, size = 1, prob = 0.5)
mcat <- as.factor(rpois(100, 1))
dummy <- as.factor(rbinom(100, size = 1, prob = 0.5))
cont <- runif(100, -10, 10)
data <- data.frame(y, mcat, dummy, cont)

# run glm logit model
glm_fit <- glm(y ~ mcat + dummy * cont, data = data, family = binomial(link = "logit"))

# compute average marginal effects
glm_ame <- ame(glm_fit, cont_vars = c("cont"), fac_vars = c("mcat", "dummy"), nsim = 1000)
```

Use the `summary()` command to display the average marginal effects and first differences of the model:
```
Continuous variables: cont 
Factor variables:     mcat, dummy 

             dydx Std. Error z value Pr(>|z|)  
cont   -0.0131441  0.0085523 -1.5369  0.12183  
mcat1  -0.3265937  0.1794777 -1.8197  0.06743 .
mcat2  -0.4169860  0.1796025 -2.3217  0.01984 *
mcat3  -0.1594161  0.1710301 -0.9321  0.34426  
mcat4  -0.3370436  0.1735345 -1.9422  0.05107 .
mcat5   0.0544205  0.1892945  0.2875  0.75826  
mcat6  -0.1293696  0.2022604 -0.6396  0.51197  
mcat7  -0.2234591  0.5102121 -0.4380  0.64818  
mcat8  -0.2191224  0.5145004 -0.4259  0.65678  
dummy1  0.0544613  0.0841782  0.6470  0.50729  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Note: For a factor variable f, dydx corresponds to the first difference E(Y|f_i) - E(Y|f_0)
```
#### 2. Linear model mtcars-data
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
