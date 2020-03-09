# (PART) Quantifying Uncertainty {-} 

# The Bootstrap {#bootstrap-main}

## Introduction

* The jackknife and the bootstrap are nonparametric procedures that are mainly used for finding standard errors
and constructing confidence intervals.

**Why use the bootstrap?**

1. To find better standard errors and/or confidence intervals when the standard approximations do not work very well.
1. To find standard errors and/or confidence intervals when you have no idea how to compute reasonable standard errors.

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

**Example: Inference for $e^{\mu}$**

* Suppose we have i.i.d. data $X_{1}, \ldots, X_{n} \sim \textrm{Logistic}( \mu, s)$,
meaning that $E(X_{i}) = \mu$ and $\textrm{Var}(X_{i}) = \sigma^{2} = s^{2}\pi^{2}/3$.

* Suppose our goal is to construct a confidence interval for the parameter $\theta = e^{\mu}$.

* The traditional approach to constructing a confidence interval uses the fact that
\begin{equation}
\sqrt{n}\Big( e^{\bar{X}} - e^{\mu} \Big) \longrightarrow \textrm{Normal}(0, \sigma^{2}e^{2\mu}) \nonumber
\end{equation}
so that we can assume $e^{\bar{X}}$ has a roughly Normal distribution with mean $e^{\mu}$ 
and standard deviation $\sigma e^{\mu}/\sqrt{n}$. This approximation is based on 
a Central Limit Theorem and "delta method" argument.

* Using this Normal approximation for $e^{\bar{X}}$, the $95\%$ confidence interval for
$e^{\mu}$ is
\begin{equation}
\Big[ e^{\bar{X}} - 1.96 \times \frac{\hat{\sigma}e^{\bar{X}}}{\sqrt{n}},
e^{\bar{X}} + 1.96 \times \frac{\hat{\sigma}e^{\bar{X}}}{\sqrt{n}}
\Big]
(\#eq:normal-approx-emu)
\end{equation}

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* For specific choices of $n$, how good is the Normal approximation for the distribution of $e^{\bar{X}}$?

* The figure below shows a histogram for the simulated distribution of $e^{\bar{X}}$ when $n=50$.
The density of the Normal approximation is also shown in this figure.

![(\#fig:unnamed-chunk-1)Histogram of simulated values of exp(sample mean) with density of the Normal approximation overlaid. This assumes n=50.](09-bootstrapLatex_files/figure-latex/unnamed-chunk-1-1.pdf) 

* As we can see from the above histogram, the Normal approximation is not terrible. However, 
it really does not capture the skewness in the distribution of $e^{\bar{X}}$ correctly.

* This could effect the coverage performance of confidence intervals which use
Normal approximation \@ref(eq:normal-approx-emu).

* The bootstrap offers an alternative approach for constructing confidence
intervals which does not depend on parametric approximations such as \@ref(eq:normal-approx-emu).

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

**Example: Inference for the Correlation**

* The sample correlation $\hat{\rho}$ which estimates the correlation $\rho = \textrm{Corr}(X_{i}, Y_{i})$
between $X_{i}$ and $Y_{i}$ is defined as
\begin{equation}
\hat{\rho} = \frac{\sum_{i=1}^{n}(X_{i} - \bar{X})(Y_{i} - \bar{Y})}{\sqrt{\sum_{i=1}^{n}(X_{i} - \bar{X})^{2}}\sqrt{\sum_{i=1}^{n}(Y_{i} - \bar{Y})}} \nonumber
\end{equation}

* Even such a relatively straightforward estimate has a pretty complicated formula for the standard error if 
you use a multivariate delta method argument:
\begin{equation}
\textrm{std.err}_{corr}
= \Bigg\{ \frac{\hat{\rho}^{2}}{4n}\Bigg[ \frac{\hat{\mu}_{40}}{\hat{\mu}_{20}^{2}} + \frac{\hat{\mu}_{04}}{\hat{\mu}_{02}^{2}} + \frac{2\hat{\mu}_{22}}{\hat{\mu}_{20}\hat{\mu}_{02} } + \frac{4\hat{\mu}_{22}}{\hat{\mu}_{11}^{2}} - \frac{4\hat{\mu}_{31}}{\hat{\mu}_{11}\hat{\mu}_{20} } + - \frac{4\hat{\mu}_{13}}{\hat{\mu}_{11}\hat{\mu}_{02} } \Bigg] \Bigg\}^{1/2}
(\#eq:rho-stderr)
\end{equation}
where
\begin{equation}
\hat{\mu}_{hk} = \sum_{i=1}^{n}(X_{i} - \bar{X})^{h}(Y_{i} - \bar{Y})^{k} \nonumber
\end{equation}

* Another popular approach for constructing a confidence interval is to use Fisher's "z-transformation"
\begin{equation}
z = \frac{1}{2} \ln\Big( \frac{1 + \hat{\rho}}{1 - \hat{\rho}}  \Big)
(\#eq:rho-stderr-ztrans)
\end{equation}
where it is argued that $z$ has a roughly Normal distribution with mean
$\tfrac{1}{2}\ln\{ (1 + \rho)/(1 - \rho) \}$ and standard deviation $1/\sqrt{n - 3}$.

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* The bootstrap allows to totally bypass the need to derive tedious formulas for the standard error such as
\@ref(eq:rho-stderr) or bypass the need to use clever transformations such as \@ref(eq:rho-stderr-ztrans). 

* For many more complicated estimates deriving formulas such as \@ref(eq:rho-stderr) or transformations such 
as \@ref(eq:rho-stderr-ztrans) may not even be feasible.

* The bootstrap provides an automatic way of constructing confidence intervals. You only 
need to be able to compute the estimate of interest.


## Description of the Bootstrap

* Suppose we have a statistic $T_{n}$ that is an estimate of some quantity of interest $\theta$.

* For example:
    + $T_{n}$ - sample mean, $\theta = E(X_{1})$.
    + $T_{n}$ - sample correlation, $\theta = \textrm{Corr}(X_{1}, Y_{1})$.
    + $T_{n}$ - sample median, $\theta = F^{-1}(1/2)$; that is, the true median. 

* Typically, $T_{n}$ can be represented as a function of our sample $X_{1}, \ldots, X_{n}$
\begin{equation}
T_{n} = h\Big( X_{1}, \ldots, X_{n}   \Big)  \nonumber
\end{equation}

* Suppose we want to estimate the standard deviation of $T_{n}$. 

* An estimate of the standard deviation of $T_{n}$ is referred to as the **standard error**. 

*Confidence intervals are often based on subtracting or adding the standard error, 
e.g.
\begin{equation}
CI = T_{n} \pm z_{\alpha/2} \times \textrm{standard error}  \nonumber 
\end{equation}

* The bootstrap estimates the standard deviation of $T_{n}$ by repeatedly subsampling from 
the original data and computing the value of the statistic $T_{n}$ on each subsample. 

* More generally, we can use the bootstrap not just to find the variance of $T_{n}$
but to characterize the distribution of $T_{n}$.

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

**The Bootstrap Procedure**

* In our description of the bootstrap, we will assume that we have the following ingredients:
    + $\mathbf{X} = (X_{1}, \ldots, X_{n})$ where $X_{1}, \ldots, X_{n}$ are i.i.d. observations.
    + The statistic $T_{n}$ of interest. $T_{n} = h(X_{1}, \ldots, X_{n})$.
    + $T_{n}$ is an estimate of $\theta$.
    + $R$ - the number of bootstrap replications.

* The bootstrap works in the following way:
* For $r = 1, \ldots, R$:
    + Draw a sample of size $n$ $(X_{1}^{*}, \ldots, X_{n}^{*})$ by sampling with replacement from $\mathbf{X}$.
    + Compute $T_{n,r}^{*} = h(X_{1}^{*}, \ldots, X_{n}^{*})$. 

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* We will refer to each sample $(X_{1}^{*}, \ldots, X_{n}^{*})$ as a **bootstrap sample**.

* We will refer to $T_{n,r}^{*}$ as a **bootstrap replication**.

* The bootstrap estimate for the standard error of $T_{n}$ is
\begin{equation}
se_{boot} = \Bigg[ \frac{1}{R-1} \sum_{r=1}^{R} \Big( T_{n,r}^{*} - \frac{1}{R} \sum_{r=1}^{R} T_{n,r}^{*} \Big)^{2} \Bigg]^{1/2} \nonumber
\end{equation}

* We can even use our bootstrap replications to get an approximation $\hat{G}_{n}^{*}(t)$ for the cumulative distribution function 
$G_{n}(t) = P(T_{n} \leq t)$ of $T_{n}$:
\begin{equation}
\hat{G}_{n}^{*}(t) = \frac{1}{R} \sum_{r=1}^{R} I\Big( T_{n,r}^{*} \leq t \Big) \nonumber
\end{equation}

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* The **normal bootstrap standard error confidence interval** is defined as 
\begin{equation}
\Big[ T_{n} - z_{\alpha/2} se_{boot}, T_{n} + z_{\alpha/2}se_{boot} \Big] \nonumber
\end{equation}

* The **bootstrap percentile confidence interval** uses the percentiles of the
boostrap replications $T_{n,1}^{*}, \ldots, T_{n,R}^{*}$ to form a confidence interval.

* The bootstrap $100 \times \alpha/2$ and $100 \times (1 - \alpha/2)$ percentiles are roughly defined as
\begin{eqnarray}
T_{[\alpha/2]}^{boot} &=& \textrm{the point } t^{*} \textrm{ such that } 100\alpha/2 \textrm{ of the bootstrap replications are less than } t^{*} \nonumber \\
T_{1 - [\alpha/2]}^{boot} &=& \textrm{the point } t^{*} \textrm{ such that } 100\alpha/2 \textrm{ of the bootstrap replications are less than } t^{*} \nonumber 
\end{eqnarray}

* The level $100 \times (1 - \alpha) \%$ level boostrap percentile confidence interval 
is then
\begin{equation}
\Big[ T_{[\alpha/2]}^{boot}, T_{[1 - \alpha/2]}^{boot} \Big]  \nonumber
\end{equation}

* More precisely, the bootstrap percentiles are obtained by looking at the inverse of the estimated cdf of $T_{n}$
\begin{equation}
T_{[\alpha/2]}^{boot} = \hat{G}_{n}^{*, -1}(\alpha/2)  \qquad  T_{[1 - \alpha/2]}^{boot} = \hat{G}_{n}^{*, -1}(1 - \alpha/2) \nonumber
\end{equation}

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* The bootstrap approach for computing standard errors and confidence intervals
is very appealing due to the fact that it is **automatic**.

* That is, we do not expend any effort deriving formulas for the variance of $T_{n}$
and/or making asymptotic arguments for the distribution of $T_{n}$.

* We only need to be able to compute $T_{n}$ many times, and the bootstrap procedure
will automatically produce a confidence interval for us.


### Example: Confidence Intervals for the Rate Parameter of an Exponential Distribution

* Suppose we have i.i.d. data $X_{1}, \ldots, X_{n}$ from an Exponential distribution with rate parameter $1/\lambda$.
That is, the pdf of $X_{i}$ is
\begin{equation}
f(x)
= \begin{cases}
\frac{1}{\lambda}e^{-x/\lambda} & \text{  if  }  x > 0 \nonumber \\
0  & \text{ otherwise }
\end{cases}
\end{equation}

* This means that 
\begin{equation}
E( X_{i} )  = \lambda \quad \textrm{ and } \quad \textrm{Var}( X_{i} ) = \lambda^{2} \nonumber
\end{equation}

* If using the usual Normal approximation for constructing a confidence interval for $\lambda$, you would
rely on the following asymptotic result:
\begin{equation}
\frac{\sqrt{n}(\bar{X} - \lambda)}{ \bar{X} }  \longrightarrow \textrm{Normal}(0, 1) \nonumber
\end{equation}

* In other words, for large $n$, $\bar{X}$ has an approximately Normal distribution with mean $\lambda$
and standard deviation $\bar{X}/\sqrt{n}$.


* The standard error in this case is $\bar{X}/\sqrt{n}$, and a $95\%$ confidence interval for $\lambda$ is 
\begin{equation}
\Bigg[ \bar{X} - 1.96 \times \frac{\bar{X}}{\sqrt{n}}, \bar{X} + 1.96 \times \frac{\bar{X}}{\sqrt{n}} \Bigg] \nonumber
\end{equation}

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* Let's do a small simulation to see how the Normal approximation confidence interval compares with
bootstrap-based confidence intervals.

* We will compare the Normal-approximation confidence interval with both the standard error bootstrap
confidence interval and the percentile bootstrap confidence interval. 




```r
xx <- rexp(50, rate=2) ## data
R <- 500   ## number of bootstrap replications
boot.mean <- rep(0, R)
for(k in 1:R) {
   boot.samp <- sample(1:50, size=50, replace=TRUE)
   xx.boot <- xx[boot.samp]   ## this is the bootstrap sample
   boot.mean[k] <- mean(xx.boot)  ## this is the kth bootstrap replication
}
```


```r
par.ci <- c(mean(xx) - 1.96*mean(xx)/sqrt(50), mean(xx) + 1.96*mean(xx)/sqrt(50))
boot.ci.sd <- c(mean(xx) - 1.96*sd(boot.mean), mean(xx) + 1.96*sd(boot.mean))
boot.ci.quant <- quantile(boot.mean, probs=c(.025, .975))
```

* The normal-approximation confidence interval is

```r
round(par.ci, 2)
```

```
## [1] 0.31 0.54
```
* The standard error boostrap confidence interval is

```r
round(boot.ci.sd, 2)
```

```
## [1] 0.32 0.52
```
* The percentile bootstrap confidence interval

```r
round(boot.ci.quant, 2)
```

```
##  2.5% 97.5% 
##  0.32  0.52
```

![](09-bootstrapLatex_files/figure-latex/unnamed-chunk-8-1.pdf)<!-- --> 

## Why is the Bootstrap Procedure Reasonable?



