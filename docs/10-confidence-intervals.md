# The Bootstrap and Confidence Intervals {#ci}
 

```r
x <- rnorm(100)
tmp <- ar(x, aic=FALSE, order.max = 1)
tmp$asy.var.coef  ## asymptotic variance of alpha.hat
```

```
##            [,1]
## [1,] 0.01001363
```


### Using the Bootstrap in Regression

* In linear regression with a single, univariate covariate, we work with the following model
\begin{equation}
Y_{i} = \beta_{0} + \beta_{1}x_{i} + \varepsilon_{i}, \qquad i = 1, \ldots, n.  \nonumber 
\end{equation}
    + $Y_{i}$ - the responses
    + $x_{i}$ - the covariates
    + $\beta_{0}, \beta_{1}$ - the regression coefficients
    + $\varepsilon_{i}$ - the residuals
    
* Typically, confidence intervals for the regression coefficients $\beta_{0}$ and $\beta_{1}$
are constructed under the assumption that $\varepsilon_{i} \sim \textrm{Normal}(0, \sigma^{2})$.


  
