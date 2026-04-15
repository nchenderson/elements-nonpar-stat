# Conformal Prediction {#conformal-prediction}

---


## Confidence Intervals and Prediction Intervals

* **Confidence Intervals** are intervals constructed to **"cover"** the parameter of a statistical model.
    + The mean $\mu$ from a parametric model.
    + Terms such as $F(1)$ from a nonparametric model where $F$ is a cdf.

* In this section, we will consider the construction of **prediction intervals**.
    + These are typically used in the context of regression models.

---

- Let's a consider a regression model for data $(Y_{i}, x_{i})$ for $i = 1, \ldots, n$ 
that has an intercept and a slope:
\begin{equation}
Y_{i} = \beta_{0} + \beta_{1}x_{i} + \varepsilon_{i},
\end{equation}
where $\varepsilon_{i} \sim \textrm{Normal}(0, \sigma^{2})$.

- The usual $95\%$ **confidence interval** for $\beta_{1}$ is:
\begin{equation}
\hat{CI}(Y_{1}, \ldots, Y_{n}) = \Big[ \hat{\beta}_{1} - 1.96 \times \frac{\hat{\sigma}}{\sqrt{\sum_{i}x_{i}^{2} - n\bar{x}^{2}}}, \hat{\beta}_{1} + 1.96 \times \frac{\hat{\sigma}}{\sqrt{\sum_{i}x_{i}^{2} - n\bar{x}^{2}}} \Big]
\end{equation}

- This has $95\%$ **"coverage"** of $\beta_{1}$ in the sense that
\begin{equation}
P\Big\{ \beta_{1} \in \hat{CI}(Y_{1}, \ldots, Y_{n}) \Big\} = 0.95
\end{equation}

- Practically, if you imagine that you had **replicated** outcomes $Y_{1s}, \ldots, Y_{ns}$ for $s = 1, \ldots, S$
generated from the same model, you should expect that your series 
of constructed confidence intervals $\hat{CI}(Y_{1s}, \ldots, Y_{ns})$ should satisfy:
\begin{equation}
\frac{1}{S}\sum_{j=1}^{S} I\Big( \beta_{1} \in \hat{CI}(Y_{1s}, \ldots, Y_{ns}) \Big) \approx 1 - \alpha
\end{equation}

---

- A quick simulation which **confirms the coverage property** of the confidence intervals for the regression coefficients is:

``` r
nsims <- 1000 ## number of simulated data sets
n <- 100
xx <- runif(n)
beta0 <- 0.5
beta1 <- 1

cover <- rep(NA, nsims)
for(s in 1:nsims) {
  Y <- beta0 + beta1*xx + rnorm(n)
  
  mod_fit <- lm(Y ~ xx)
  ## Use method confint to get 95% confidence intervals for regression coefficients:
  lower_ci <- confint(mod_fit)[2,1]
  upper_ci <- confint(mod_fit)[2,2]  
    
  cover[s] <- ifelse(beta1 >= lower_ci & beta1 <= upper_ci, 1, 0)
}
mean(cover) ## Should be close to 0.95
```

```
## [1] 0.954
```

---

* **Prediction intervals** for regression have a different interpretation than confidence intervals.

* For prediction intervals, you usually imagine a **"future observation"** $Y_{n+1}$ that come from 
the same regression model
\begin{equation}
Y_{n+1} = \beta_{0} + \beta_{1}x_{n+1} + \varepsilon_{n+1},
\end{equation}


* A $100 \times (1 - \alpha)$ **prediction interval** $\hat{C}(x_{n+1})$ constructed from data $(Y_{1}, \mathbf{x}_{1}), \ldots, (Y_{n}, \mathbf{x}_{n})$
is supposed to have the following property:
\begin{equation}
P\Big( Y_{n+1} \in \hat{C}(x_{n+1}) \Big) = 1 - \alpha
\end{equation}

* For the linear regression model with a single covariate, the standard $95\%$ prediction interval has the form
\begin{equation}
\hat{\beta}_{0} + \hat{\beta}_{1}x_{n+1} \pm 1.96 \times \hat{\sigma}\sqrt{1 + (1,x_{n+1})^{T}(\mathbf{X}^{T}\mathbf{X})^{-1}(1, x_{n+1})}
\end{equation}

* The above prediction interval works for linear regression.
    + We want a general procedure that works for any procedure you would use for modeling the conditional mean $E\{ Y_{i} | x_{i} \}$.


## Conformal Inference Procedure for Prediction Intervals

* The **split-sample conformal inference** procedure works for the setting where you're thinking of
outcomes $Y_{i}$ coming from the following model:
\begin{equation}
Y_{i} = f(\mathbf{x}_{i}) + \varepsilon_{i}
\end{equation}

* You **do not** need to know the "true form" of $f$.
    - You could fit a linear model, but the true $f$ could be a nonlinear function of covariates.
    - The conformal inference procedure will still work.

* $f$ can just be the fitted values returned by a **"black box" machine learning procedure**.
    - For example, $f(\mathbf{x}_{i})$ could be the fitted values returned by running boosting or random forest.

---

* Let $\mathcal{D}$ denote your dataset.
    + This contains pairs of the form $(\mathbf{x}_{i}, Y_{i})$
    + $\mathcal{D} = \{(\mathbf{x}_{1}, Y_{1}), \ldots (\mathbf{x}_{n}, Y_{n}) \}$.

1. The first step is to **split** your data $\mathcal{D}$ further into **non-overlapping** sets:
     + $\mathcal{A}_{1}$ - the indeces of a **"proper"** dataset with $n_{1}$ observations.
     + $\mathcal{A}_{2}$ - the indeces of a **"calibration"** set with $n_{2}$ observations.
     + $\mathcal{D}_{1}$ - the proper dataset: $\mathcal{D}_{1} = \{(\mathbf{x}_{i}, Y_{i}): i \in \mathcal{A}_{1} \}$.
     + $\mathcal{D}_{2}$ - the calibration dataset: $\mathcal{D}_{2} = \{(\mathbf{x}_{i}, Y_{i}): i \in \mathcal{A}_{2} \}$.
     
2. Using **only data** from $\mathcal{D}_{1}$ apply your regression/machine learning procedure to build
a function $\hat{f}_{\mathcal{D}_{1}}(\mathbf{x}_{i})$ that predicts $Y_{i}$ from $\mathbf{x}_{i}$.

3. For $i \in \mathcal{A}_{2}$, define the calibration set **absolute residuals** $R_{i}$ as
\begin{equation}
R_{i} = | Y_{i} - \hat{f}_{\mathcal{D}_{1}}(\mathbf{x}_{i})|
\end{equation}


4. Compute the $100 \times (1 - \alpha)$ quantile (usually $\alpha = 0.05$) of these calibration residuals
\begin{equation}
\hat{q}_{\mathcal{D}_{2}, \alpha} = \textrm{$1 - \alpha$ quantile from residuals} R_{i} \textrm{ such that } i \in \mathcal{A}_{2}.
\end{equation}

5. Then, use the quantile $\hat{q}_{\mathcal{D}_{2}}$ to form the conformal $(1 - \alpha)$ prediction set
\begin{equation}
\hat{C}_{n}( \mathbf{x} ) = \Big[\hat{f}_{\mathcal{D}_{1}}(\mathbf{x}) - \hat{q}_{\mathcal{D}_{2}, \alpha},\hat{f}_{\mathcal{D}_{1}}(\mathbf{x}) + \hat{q}_{\mathcal{D}_{2}, \alpha}\Big]
\end{equation}

---

* Applying the **above 5 steps** results in **conformal prediction intervals** that satisfies the following property
\begin{equation}
1 - \alpha \leq P\Big( Y_{n+1} \in \hat{C}_{n}(\mathbf{x}_{n+1}) \Big| \textrm{data in proper set}) < 1 - \alpha + \frac{1}{n_{2} + 1}
\end{equation}

* Remarkably, this interval is "**distribution free**".
      - It does not assume that one has correctly specified the model for the outcomes.

* The main assumption required are: 
     - $(Y_{n+1}, \mathbf{x}_{n+1})$ is independent from $((Y_{1}, \mathbf{x}_{1}), \ldots, (Y_{n+1}, \mathbf{x}_{n}))$
     - The joint distribution of $(Y_{n+1}, \mathbf{x}_{n+1})$ is the same as $(Y_{i}, \mathbf{x}_{i})$.

---

* As an example, let's try to generate a conformal prediction interval for a **linear regression** example:

* First generate the data

``` r
n <- 100
xx <- runif(n)
beta0 <- 0.5
beta1 <- 1

Y <- beta0 + beta1*sqrt(xx) + rnorm(n)
```
    - Note that the data are generated from the model $Y_{i} = \beta_{0} + \beta_{1}\sqrt{x_{i}} + \varepsilon_{i}$.
    
    - We will fit the **"misspecified"** model $Y_{i} = \beta_{0} + \beta_{1}x_{i} + \varepsilon_{i}$, but the conformal inference procedure should still work. 
    
    
* Split this data set into a "proper" training set and calibration set

``` r
A1 <- sample(1:n, size=50)
A2 <- setdiff(1:n, A1)

proper_dat <- data.frame(Y=Y[A1], xx=xx[A1])
calibration_dat <- data.frame(Y=Y[A2], xx=xx[A2])
```

* Using `D1`, fit a linear regression model

``` r
proper_mod <- lm(Y ~ xx, data=proper_dat)
```

* Get absolute residuals on calibration dataset

``` r
calibration_fitted <- predict(proper_mod, newdat=calibration_dat)
calibration_resids <- abs(calibration_dat$Y - calibration_fitted)
```

* Get 95th quantile of these residuals

``` r
qhat <- quantile(calibration_resids, probs=0.95)
print(qhat)
```

```
##      95% 
## 2.270708
```

* We can now use `qhat` to get prediction intervals for a "new dataset"

``` r
## Generate new dataset
n <- 100
xx_new <- runif(n)
beta0 <- 0.5
beta1 <- 1

Y_new <- beta0 + beta1*sqrt(xx) + rnorm(n)
newdataset <- data.frame(Y=Y_new, xx=xx_new)

### Construct prediction intervals as a n x 2 matrix
ConformalInterval <- matrix(NA, nrow=n, ncol=2)
ConformalInterval[,1] <- predict(proper_mod, newdat=newdataset) - qhat
ConformalInterval[,2] <- predict(proper_mod, newdat=newdataset) + qhat

print(head(ConformalInterval))
```

```
##           [,1]     [,2]
## [1,] -1.276252 3.265163
## [2,] -1.192513 3.348902
## [3,] -1.375029 3.166386
## [4,] -1.050272 3.491143
## [5,] -1.403428 3.137987
## [6,] -1.020907 3.520508
```

* Plot fitted values and prediction intervals:
<img src="11-conformalprediction_files/figure-html/unnamed-chunk-8-1.png" alt="" width="672" />

* You can check the **prediction coverage** of these intervals with the following code:

``` r
## This shouldn't be that far off 0.95, but there will be 
## considerable variability since this is not a very large dataset
mean(Y_new > ConformalInterval[,1] & Y_new < ConformalInterval[,2])
```

```
## [1] 0.98
```

## Why does this work?

* The main justification for the validity of this procedure comes from 
looking at the **calibration residuals** $R_{i}$, for $i \in \mathcal{A}_{2}$ 
and the **test residual** $R_{n + 1} = | Y_{n+1} - \hat{f}_{\mathcal{D}_{1}}(\mathbf{x}_{n+1})|$.

* Specifically, $(R_{n+1}, R_{i}, i \in \mathcal{A}_{2})$ is a collection of **i.i.d random variables** that are independent from observations in $\mathcal{D}_{1}$.
    - This is true because $\hat{f}_{\mathcal{D}_{1}}(\mathbf{x})$ was built from the **proper training set** and ...
    
    - The values of $R_{i}$, for $i \in \mathcal{A}_{2}$ only uses outcomes from the **calibration dataset**.
    

* Because of the pairs $(\mathbf{x}_{i}, Y_{i})$ are i.i.d., the probability that $R_{n+1}$ is **less than the $100(1 - \alpha)$ quantile** of the residuals is very close to $1 - \alpha$. 

* To be more specific:
\begin{eqnarray}
P\Big( Y_{n+1} \in \hat{C}_{n}(\mathbf{x}_{n+1}) \Big| \mathcal{D}_{1} \Big)
&=& P(\hat{f}_{\mathcal{D}_{1}}(\mathbf{x}_{n+1}) - \hat{q}_{\mathcal{D}_{2}, \alpha} \leq Y_{n+1} \leq \hat{f}_{\mathcal{D}_{1}}(\mathbf{x}_{n+1}) + \hat{q}_{\mathcal{D}_{2}, \alpha} \Big| \mathcal{D}_{1} \Big) \nonumber \\
&=& P(- \hat{q}_{\mathcal{D}_{2}, \alpha} \leq Y_{n+1} - \hat{f}_{\mathcal{D}_{1}}(\mathbf{x}_{n+1}) \leq \hat{q}_{\mathcal{D}_{2}, \alpha} \Big| \mathcal{D}_{1} \Big) \nonumber \\
&=& P( R_{n+1} \leq \hat{q}_{\mathcal{D}_{2}, \alpha} \Big| \mathcal{D}_{1} \Big) \nonumber \\
&\approx& 1 - \alpha
\end{eqnarray}
    
    
## Conformal Prediction Intervals with Tree Boosting.

* Boosting with trees is a method for estimating $f(\mathbf{x})$ in the following
regression model
\begin{equation}
Y_{i} = f(\mathbf{x}_{i}) + \varepsilon_{i}
\end{equation}

* You **do not** need to know the details of how $f(\mathbf{x})$ is estimated
in order to create conformal prediction intervals.

* We will just use the `gbm` package in **R** to estimate $f(\mathbf{x})$.
    + We won't discuss the details of how `gbm` estimates $f(\mathbf{x})$ now.
    + You can think of it as a "black box" procedure that generates fitted values.

---

* To illustrate conformal inference with boosting, we will first generate an example dataset with 20 covariates

``` r
n <- 2000
p <- 20
beta0 <- c(1, -1, 2, -2, rep(0, 16))

## Generate design matrix and outcomes:
X <- matrix(rnorm(n*p), nrow=n, ncol=p)
Y <- X%*%beta0 + rnorm(n)
```


* First, split this dataset into a "proper" set (of size 1000) and a calibration set

``` r
A1 <- sample(1:n, size=1000)
A2 <- setdiff(1:n, A1)

proper_dat <- data.frame(Y=Y[A1], X[A1,])
calibration_dat <- data.frame(Y=Y[A2], X[A2,])
```

* Using the proper set, use boosting to estimate f 

``` r
library(gbm)
```

```
## Loaded gbm 2.2.3
```

```
## This version of gbm is no longer under development. Consider transitioning to gbm3, https://github.com/gbm-developers/gbm3
```

``` r
gbm_mod <- gbm(Y ~ ., data = proper_dat, 
           distribution = "gaussian", n.trees = 200, cv.folds=5)

## Find the best number of trees using cross-validation
best.iter <- gbm.perf(gbm_mod, method = "cv")
```

<img src="11-conformalprediction_files/figure-html/unnamed-chunk-12-1.png" alt="" width="672" />

``` r
print(best.iter)
```

```
## [1] 198
```

``` r
## Use boosting with best number of trees
gbm_mod_final <- gbm(Y ~ ., data = proper_dat, 
                     distribution = "gaussian", n.trees = best.iter)
```

* Get **absolute residuals** on calibration dataset

``` r
calibration_fitted <- predict(gbm_mod_final, newdat=calibration_dat)
```

```
## Using 198 trees...
```

``` r
calibration_resids <- abs(calibration_dat$Y - calibration_fitted)
```

* Get **95th quantile** of these residuals in the calibration set

``` r
qhat <- quantile(calibration_resids, probs=0.95)
```

* We can now use `qhat` to get prediction intervals for a "new dataset"
    + We generated the new dataset with the same regression parameters used to generate the first dataset.

``` r
X_new <- matrix(rnorm(n*p), nrow=n, ncol=p)

beta0 <- c(1, -1, 2, -2, rep(0, 16))
Y_new <- X_new%*%beta0 + rnorm(n)

##
newdataset <- data.frame(Y=Y_new, X_new)

### Construct prediction intervals as a n x 2 matrix
ConformalInterval <- matrix(NA, nrow=n, ncol=2)
ConformalInterval[,1] <- predict(gbm_mod_final, newdat=newdataset) - qhat
```

```
## Using 198 trees...
```

``` r
ConformalInterval[,2] <- predict(gbm_mod_final, newdat=newdataset) + qhat
```

```
## Using 198 trees...
```



* Check the **prediction coverage**:

``` r
## This shouldn't be that far off 0.95, but there will be 
## considerable variability since this is not a very large dataset
mean(Y_new > ConformalInterval[,1] & Y_new < ConformalInterval[,2])
```

```
## [1] 0.948
```



---
