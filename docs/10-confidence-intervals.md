# Bootstrap Examples and the Jackknife {#ci}
 

## The Parametric Bootstrap for an AR(1) model

* Consider the time series $X_{1}, X_{2}, \ldots, X_{m}$. Here,
$X_{t}$ denotes an observation made at time $t$.

* An autoregressive model of order 1 (usually called an AR(1) model) for this time series is
\begin{eqnarray}
X_{1} &=& \frac{c_{0}}{1 - \alpha} + \varepsilon_{1} \nonumber \\
X_{t} &=& c_{0} + \alpha X_{t-1} + \varepsilon_{t}, \qquad t=2,\ldots,m. \nonumber
\end{eqnarray}

* It is usually assumed that $|\alpha| < 1$.

* In the AR(1) model, it is assumed that 
    + $E(\varepsilon_{t}) = 0$
    + $\textrm{Var}(\varepsilon_{t}) = \sigma^{2}$,
    + $\varepsilon_{2}, \ldots, \varepsilon_{m}$ are i.i.d.
    + $\varepsilon_{t}$ and $X_{t-1}$ are independent.

* In addition to these assumptions, we will assume that
\begin{equation}
\varepsilon_{t} \sim \textrm{Normal}(0, \sigma^{2})  \nonumber 
\end{equation}

* The AR(1) model implies that 
\begin{equation}
\textrm{Corr}(X_{t}, X_{t-1}) = \alpha  \nonumber 
\end{equation}
and, more generally, that
\begin{equation}
\textrm{Corr}(X_{t}, X_{t-p}) = \alpha^{p}  \nonumber
\end{equation}


---

* For known values of $c_{0}, \alpha$, and $\sigma^{2}$, we can simulate
an AR(1) time series with the following `R` code:

```r
SimulateParAR1 <- function(m, c0, alpha, sig.sq) {
     xx <- numeric(m)
     xx[1] <- c0/(1 - alpha) + rnorm(1, sd=sqrt(sig.sq))
     for(t in 2:m) { 
         xx[t] <- c0 + alpha*xx[t-1] + rnorm(1, sd=sqrt(sig.sq))
     }
     return(xx)
}
```


<img src="10-confidence-intervals_files/figure-html/unnamed-chunk-2-1.png" width="672" />
 
* In `R`, estimates of $c_{0}, \alpha,$ and $\sigma^{2}$ can be found by using the `ar` function. For example,

```r
x <- SimulateParAR1(1000, 1, 0.8, sig.sq=.25)
ar1.fit <- ar(x, aic=FALSE, order.max = 1, method="mle")

c0.est <- ar1.fit$x.mean*(1 - ar1.fit$ar)
alpha.est <- ar1.fit$ar
sigsq.est <- ar1.fit$var.pred
```

---

* Suppose we want to construct confidence intervals for $\alpha$ and $\sigma$ using a bootstrap method.

* Using the direct, nonparametric bootstrap described in the previous chapter will not work
because our observations are not independent. There are "block bootstraps" that
are designed to work for time series, but we will not discuss those here (see e.g., @buhlmann2002 or Chapter 8 of @davison1997 for 
more details).

* With the parametric bootstrap, we only have to use the following steps to generate bootstrap replications
$\hat{\alpha}_{r}^{*}$ and $\hat{\sigma}_{r}^{2,*}$ for estimates of $\alpha$ and $\hat{\sigma}^{2}$.
 
* For $r = 1, \ldots, R$:
    + Simulate a time series $X_{1}^{*}, \ldots, X_{m}^{*}$ from an AR(1) model with parameters $(\hat{c}_{0}, \hat{\alpha}, \hat{\sigma}^{2})$.
    + Compute $\hat{\alpha}_{r}^{*} = \hat{\alpha}(X_{1}^{*}, \ldots, X_{m}^{*})$.
    + Compute $\hat{\sigma}_{r}^{2,*} = \hat{\sigma}^{2}(X_{1}^{*}, \ldots, X_{m}^{*})$

---

* To see how this parametric bootstrap works, we will use the `nhtemp` dataset that is available in `R`.

<img src="10-confidence-intervals_files/figure-html/unnamed-chunk-4-1.png" width="672" />

* The `nhtemp` dataset contains the mean annual temperature in New Haven, Connecticut from the years 1912-1971

```r
head(nhtemp)
```

```
## [1] 49.9 52.3 49.4 51.1 49.4 47.9
```

* The estimated autocorrelation parameter $\alpha$ is about ? for this data

```r
ar1.temp <- ar(nhtemp, aic=FALSE, order.max = 1)
c0.hat <- ar1.temp$x.mean*(1 - ar1.temp$ar)
alpha.hat <- ar1.temp$ar
sigsq.hat <- ar1.temp$var.pred
alpha.hat
```

```
## [1] 0.3148269
```

* Now, that we have estimated all the parameter of the AR(1) model, we can run our parametric bootstrap for $\hat{\alpha}$ and $\hat{\sigma}$:

```r
R <- 500
alpha.boot <- numeric(R)
sigsq.boot <- numeric(R)
for(r in 1:R) {
  x <- SimulateParAR1(60, c0=c0.hat, alpha=alpha.hat, sig.sq=sigsq.hat)
  ar1.fit <- ar(x, aic=FALSE, order.max = 1)
  
  alpha.boot[r] <- ar1.fit$ar
  sigsq.boot[r] <- ar1.fit$var.pred
}
```

* Normal bootstrap standard error confidence intervals for $\alpha$ and $\sigma^{2}$ are

```r
round(c(alpha.hat - 1.96*sd(alpha.boot), alpha.hat + 1.96*sd(alpha.boot)), 3)
```

```
## [1] 0.073 0.557
```

```r
round(c(sigsq.hat - 1.96*sd(sigsq.boot), sigsq.hat + 1.96*sd(sigsq.boot)), 3)
```

```
## [1] 0.933 2.003
```

* We can compare our confidence interval for $\alpha$ with the confidence interval
obtained from using a large-sample approximation:

```r
asymp.se <- sqrt(ar1.temp$asy.var.coef)
round(c(alpha.hat - 1.96*asymp.se, alpha.hat + 1.96*asymp.se), 3)
```

```
## [1] 0.071 0.559
```


## Using the Bootstrap in Regression

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

* The bootstrap allows us to compute confidence intervals for $(\beta_{0}, \beta_{1})$ without
relying on this normality assumption.

* How to compute bootstrap confidence intervals for $\beta_{0}$ and $\beta_{1}$?

---

* The least-squares estimates of $\beta_{0}$ and $\beta_{1}$ are
\begin{equation}
\hat{\beta}_{0} = \bar{y} - \hat{\beta}_{1}\bar{x} \qquad \qquad \hat{\beta}_{1} = \frac{\sum_{i=1}^{n}(x_{i} - \bar{x})(y_{i} - \bar{y})}{S_{xx}}  \nonumber
\end{equation}
where $S_{xx} = \sum_{i=1}^{n}( x_{i} - \bar{x})^{2}$.

* Assuming the covariates are fixed design points, the variance of $\hat{\beta}_{0}$ and $\hat{\beta}_{1}$ are
\begin{equation}
\textrm{Var}(\hat{\beta}_{0}) = \sigma^{2}\Big(\frac{1}{n} + \frac{\bar{x}}{S_{xx}} \Big) \qquad \textrm{Var}(\hat{\beta}_{1}) = \frac{\sigma^{2}}{S_{xx}} \nonumber
\end{equation}


---

* With a parametric bootstrap, we simulate outcomes $Y_{i}$ from the model
\begin{equation}
Y_{i} = \hat{\beta}_{0} + \hat{\beta}_{1}x_{i} + \varepsilon_{i}, \qquad \textrm{Normal}(0, \hat{\sigma}^{2}) \nonumber
\end{equation}
where $\hat{\beta}_{0}$ and $\hat{\beta}_{1}$ are the least-squares estimates 
and $\hat{\sigma}^{2} = \tfrac{1}{n-2}\sum_{i=1}^{n} (Y_{i} - \hat{\beta}_{0} - \hat{\beta}_{1})^{2}$. 




```r
kidney <- read.table("https://web.stanford.edu/~hastie/CASI_files/DATA/kidney.txt", 
                     header=TRUE)
```

<img src="10-confidence-intervals_files/figure-html/unnamed-chunk-11-1.png" width="672" />




  
