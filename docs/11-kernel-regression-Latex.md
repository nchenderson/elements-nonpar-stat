# (PART) Nonparametric Regression: Part I {-}

# Kernel Regression and Local Regression
 
## Introduction

* In regression we are interested in characterizing, in some way, the relationship
between a collection of responses $Y_{1},\ldots,Y_{n}$ and covariate vectors
$(\mathbf{x}_{1}, \ldots, \mathbf{x}_{n})$.

* Linear regression is one way of approaching this problem. This assumes
the expectation of $Y_{i}$ can be expressed as a linear combination
of the covariates:
\begin{equation}
E(Y_{i}| \mathbf{x}_{i}) = \beta_{0} + \sum_{j=1}^{p} x_{ij}\beta_{j}  \nonumber
\end{equation}

* More generally, we can consider the following model
\begin{equation}
Y_{i} = m( \mathbf{x}_{i} ) + \varepsilon_{i}  \nonumber
\end{equation}
    + $m(\mathbf{x}_{i})$ - the "mean function" or "regression function"
    + $\mathbf{x}_{i} = (x_{i1}, \ldots, x_{ip})$ - the $i^{th}$ covariate vector
    
* The residuals $\varepsilon_{1}, \ldots, \varepsilon_{n}$ are assumed to 
be i.i.d. and have mean zero. 

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* In a nonparametric approach, we will try to estimate $m(\mathbf{x})$ without
making any strong assumptions about the form of $m( \mathbf{x} )$. 

* The regression function $m(\mathbf{x})$ can be thought of as the function which returns 
the expectation of $Y_{i}$ given that $\mathbf{x}_{i} = \mathbf{x}$
\begin{equation}
m(\mathbf{x} ) = E(Y_{i}|\mathbf{x}_{i}=\mathbf{x}) 
\end{equation}

* Let
    + $f_{Y|X}(y|\mathbf{x})$ denote the conditional density of $Y_{i}$ given $\mathbf{x}_{i}$.
    + $f_{Y,X}(y, \mathbf{x})$ denote the joint density of $(Y_{i}, \mathbf{x}_{i})$
    + $f_{X}(\mathbf{x})$ denote the density of $\mathbf{x}_{i}$
    
* We can express the regression function as
\begin{equation}
m(\mathbf{x}) = \int_{-\infty}^{\infty} y f_{Y|X}(y|\mathbf{x}) dy = \frac{\int y f_{Y,X}(y, \mathbf{x}) dy}{ f_{X}(\mathbf{x})  }  \nonumber
\end{equation}

## Kernel Regression

* In this section, we will assume that the covariates are univariate. 
That is, $p=1$ and $\mathbf{x}_{i} = x_{i}$ where $x_{i}$ is a real number.

### The Regressogram

* The regressogram is an estimate of the mean function $m(x)$ which is 
has many similarities in its construction to the histogram.

* Similar to how we constructed the histogram, let us think about an estimate $m(x)$
that will be constant within each of a series of bins $B_{1}, \ldots, B_{D_{n}}$
\begin{eqnarray}
B_{1} &=& [ x_{0}, x_{0} + h_{n})  \nonumber \\
B_{2} &=& [x_{0} + h_{n}, x_{0} + 2h_{n})  \nonumber \\
&\vdots&  \nonumber \\
B_{D_{n}} &=& [x_{0} + (D_{n} - 1)h_{n}, x_{0} + D_{n}h_{n})  \nonumber
\end{eqnarray}

* Suppose we want to estimate $m(x)$, where $x$ belongs to the $k^{th}$ bin.
A direct estimate of this is the average of the $Y_{i}'s$ among those
$x_{i}'s$ which fall into the $k^{th}$ bin.

* Specifically, if $x \in B_{k}$, then we estimate $m(x)$ with
\begin{equation}
\hat{m}_{h_{n}}^{R}(x) =  \frac{ \sum_{i=1}^{n} Y_{i} I\big( x_{i} \in B_{k} \big) }{ \sum_{i=1}^{n} I\big( x_{i} \in B_{k} \big) } 
= \frac{1}{n_{k,h_{n}}} \sum_{i=1}^{n} Y_{i} I\big( x_{i} \in B_{k} \big), \nonumber
\end{equation}
where $n_{k,h_{n}}$ is the number of $x_{i}$ that fall into the $k^{th}$ bin when using bin width $h_{n}$.

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* The estimate $\hat{m}_{h_{n}}^{R}(x)$ of the regression function is called the **regressogram**.

* The intuition for this estimate is: if $x \in B_{k}$,
then taking an average of the reponses for $x_{i}$ in a small bin containing $x$ 
should give us a reasonable approximation for the expectation of $Y_{i}$ given that $x_{i} = x$.

* Another way of looking at the regressogram is to note that if we think of the $x_{i}$ as random variables, then for $x \in B_{k}$
\begin{eqnarray}
E\Big\{ \frac{1}{n} \sum_{i=1}^{n} Y_{i} I\big( x_{i} \in B_{k} \big) \Big\}
&=& E\Big\{  Y_{1} I\big( x_{1} \in B_{k} \big) \Big\}  \nonumber \\
&=& \int_{-\infty}^{\infty} \int_{x_{0} + (k-1)h_{n}}^{x_{0} + kh_{n}} y f_{Y,X}(y, t) dt dy  \nonumber \\
&\approx& h_{n} \int_{-\infty}^{\infty} y f_{Y,X}(y, x) dy
(\#eq:regressogram-numerator)
\end{eqnarray}
and, similarly, 
\begin{eqnarray}
E\Big\{ \frac{1}{n} \sum_{i=1}^{n}  I\big( x_{i} \in B_{k} \big) \Big\}
&=& E\Big\{  I\big( x_{1} \in B_{k} \big) \Big\}  \nonumber \\
&=& \int_{x_{0} + (k-1)h_{n}}^{x_{0} + kh_{n}}  f_{X}(t) dt  \nonumber \\
&\approx& h_{n} f_{X}(x) 
(\#eq:regressogram-denominator)
\end{eqnarray}

* Equations \@ref(eq:regressogram-numerator) and \@ref(eq:regressogram-denominator) suggest that $\hat{m}_{h_{n}}^{R}(x)$
should be a reasonable estimate of the ratio 
\begin{equation}
\int_{-\infty}^{\infty} y f_{Y,X}(y, x) dy \big/ f_{X}(x) \nonumber 
\end{equation}

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

![(\#fig:unnamed-chunk-1)Framingham Data. Regressogram estimate for a regression model with diastolic blood pressure as the response and age as the covariate. Ages from 31-71 were separated into bins of width 5 years.](11-kernel-regression-Latex_files/figure-latex/unnamed-chunk-1-1.pdf) 

![](11-kernel-regression-Latex_files/figure-latex/unnamed-chunk-2-1.pdf)<!-- --> 

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* **Exercise 11.1** Let 
\begin{equation}
\hat{\mathbf{m}} = \big( \hat{m}_{h_{n}}^{R}(x_{1}), \ldots, \hat{m}_{h_{n}}^{R}(x_{n}) \big)  \nonumber
\end{equation}
denote the vector of "fitted values" from a regressogram estimate that has $D_{n}$ bins. 
If $\mathbf{Y} = (Y_{1}, \ldots, Y_{n})$, show that you can express $\hat{\mathbf{m}}$ as 
\begin{equation}
\hat{\mathbf{m}} = \mathbf{A}\mathbf{Y}, \nonumber
\end{equation}
for an appropriately chosen $n \times n$ matrix $\mathbf{A}$.
What is the value of $\textrm{tr}(\mathbf{A})$?

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}


### The Local Average Estimator

* The regressogram can be thought of as a regression analogue of the histogram.

* The local average estimator can be thought of as a regression analogue of the
"box-type" density estimator that we described in Chapter 8.

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* For each point $x$, we are going to use a regression function estimate which
has a bin "centered" at $x$.

* Specifically, for each $x$, we will form a bin of width $2h_{n}$ around
$x$ and compute the mean of the $Y_{i}$ among those observations where the $x_{i}$ fall
into this bin. 

* In other words, we are computing an average of the $Y_{i}$ in a small region
around $x$. 

* The local average estimator $\hat{m}_{h_{n}}^{loc}(x)$ at $x$ is defined as:
\begin{eqnarray}
\hat{m}_{h_{n}}^{loc}(x) &=&
\frac{ \sum_{i=1}^{n} Y_{i}I\big( x - h_{n} < x_{i} < x + h_{n} \big) }{ \sum_{i=1}^{n} I\big( x - h_{n} < x_{i} < x + h_{n} \big) } \nonumber \\
&=& \frac{1}{n_{h_{n}}(x)} \sum_{i=1}^{n} Y_{i}I\big( x - h_{n} < x_{i} < x + h_{n} \big) \nonumber 
\end{eqnarray}
where $n_{h_{n}}(x) = \sum_{i}^{n} I\big( x - h_{n} < x_{i} < x + h_{n} \big)$.

* The local average estimator does not need to have a 
constant value within each of a few pre-specified bins.

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* We can also express the local average estimator in the following way:
\begin{equation}
\hat{m}_{h_{n}}^{loc}(x)
= \frac{\sum_{i=1}^{n} Y_{i} w\Big( \frac{x - X_{i}}{h_{n}} \Big)}{\sum_{i=1}^{n} w\Big( \frac{x - X_{i}}{h_{n}} \Big)}, 
(\#eq:local-average-estimator)
\end{equation}
where $w(t)$ is the "box" function defined as
\begin{equation}
w(t) = 
\begin{cases}
\frac{1}{2} & \textrm{ if } |t| < 1 \nonumber \\
0 & \textrm{ otherwise}  \nonumber
\end{cases}
\end{equation}

* While a local average estimate will not be a "step function" like the regressogram, the local average
estimate will typically be non-smooth and have a jagged appearance.

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* Like kernel density estimation, there is a bias/variance tradeoff to the choice of $h_{n}$.

* Smaller values of $h_{n}$ usually imply higher variance because you will be taking an
average over a relatively small number of observations.

* Larger values of $h_{n}$ usually imply higher bias because 
you will be esitmating $m(x)$ by averaging over a wide range
of $x_{i}$ values, and $m(x)$ could vary substantially over this range of $x_{i}$ values.

* Our experience in Chapter 8 suggests that we can get a smoother estimate of 
the regression if we simply replace the "box function" $w(t)$ in \@ref(eq:local-average-estimator) with a
smoother kernel function $K(t)$.

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* `R` code for computing a local average estimate $\hat{m}_{2}^{loc}(x)$ at the 
points $x = 31, 32, 33, ...., 71$ is given below

```r
xseq <- seq(31, 71, by=1)
hn <- 2
nx <- length(xseq)
m.hat.loc <- numeric(nx)
for(k in 1:nx) {
    in.bin <- framingham$age > xseq[k] - hn & framingham$age < xseq[k] + hn
    m.hat.loc[k] <- mean(framingham$diaBP[in.bin])
}

plot(framingham$age, framingham$diaBP, las=1, ylab="Diastolic Blood Pressure", 
     xlab="Age", main="Local Average Estimate with hn=2", type="n")
points(framingham$age, framingham$diaBP, pch=16, cex=0.7)
lines(xseq, m.hat.loc, lwd=3, col="red")
```

![](11-kernel-regression-Latex_files/figure-latex/unnamed-chunk-3-1.pdf)<!-- --> 

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* Let's also look at a local average estimate of the regression function for the bone mineral density dataset.

* The responses in this dataset are relative changes in the bone mineral density of adolescents.

* Specifically, reponses $Y_{i}$ and covariates $x_{i}$ are defined as
\begin{eqnarray}
Y_{i} &=& \frac{\textrm{Mineral Density at Visit 2}_{i} - \textrm{Mineral Density at Visit 1}_{i}}{\tfrac{1}{2}(\textrm{Mineral Density at Visit 2}_{i} + \textrm{Mineral Density at Visit 1}_{i})}  \nonumber \\
x_{i} &=& \frac{1}{2}(\textrm{Age at Visit 2}_{i} + \textrm{Age at Visit 1}_{i})  \nonumber
\end{eqnarray}


```r
tmp <- read.table("https://web.stanford.edu/~hastie/ElemStatLearn/datasets/bone.data", 
                  header=TRUE)
bonedat <- tmp[!duplicated(tmp$idnum),]  ## only keep the first observation of a person
```


```r
xseq <- seq(9.4, 25.2, by=.1)
hn <- 1
nx <- length(xseq)
m.hat.loc <- numeric(nx)
for(k in 1:nx) {
    in.bin <- bonedat$age > xseq[k] - hn & bonedat$age < xseq[k] + hn
    m.hat.loc[k] <- mean(bonedat$spnbmd[in.bin])
}

plot(bonedat$age, bonedat$spnbmd, las=1, ylab="Relative Change in Bone MD", 
     xlab="Age", main="Bone Data: Local Average Estimate with hn=1", type="n")
points(bonedat$age, bonedat$spnbmd, pch=16, cex=0.7)
lines(xseq, m.hat.loc, lwd=3, col="red")
abline(0, 0)
```

![](11-kernel-regression-Latex_files/figure-latex/unnamed-chunk-5-1.pdf)<!-- --> 

![](11-kernel-regression-Latex_files/figure-latex/unnamed-chunk-6-1.pdf)<!-- --> 

### k-Nearest Neighbor (k-NN) Regression

* k-nearest neighbor regression is fairly similar to the local average estimator
of the regression function.

* With k-NN, we still estimate the regression function at a particular point
by taking a type of local average around this point.

* However, k-NN takes the average over the k "nearest observations" to $x$
rather than taking an average over all the observations which fall into
a bin centered at $x$.

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* The k-NN estimator of the regression function $\hat{m}_{k}^{kNN}(x)$ is defined as
\begin{equation}
\hat{m}_{k}^{kNN}(x) = \frac{1}{k}\sum_{i=1}^{n} y_{i} I\big( x \in N_{k}(x) \big) \nonumber
\end{equation}

* Here, $N_{k}(x)$ is defined as the set of the k $x_{i}'s$ which are
closest to $x$. That is, $N_{k}(x)$ is the set of the k "nearest neighbors"
to $x$. 

* Mathematically, if we define 
\begin{equation}
d_{i}(x) = |x_{i} - x| \nonumber
\end{equation}
and order them so that $d_{(1)}(x) \leq d_{(2)}(x) \leq \ldots \leq d_{(n)}(x)$.
Then, the k nearest neighbors of $x$ would be those 
observations which correspond to the $d_{(1)}(x)$ through $d_{(k)}(x)$.

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* Like the local average estimator, increasing the value of $k$
will increase the bias of the k-NN regression function estimate while 
decreasing the value of $k$ will increase the variance of the k-NN regression function estimate.

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* **Exercise 8.2** Suppose $n=6$ and that we have the following covariate values and responses
\begin{eqnarray}
(x_{1}, x_{2}, x_{3}, x_{4}, x_{5}, x_{6}) &=& (1/7, 2/7, 3/7, 4/7, 5/7, 6/7)  \nonumber \\
(Y_{1}, Y_{2}, Y_{3}, Y_{4}, Y_{5}, Y_{6}) &=& (1.4, 0.7, 1.1, 1.3, 0.9, 1.7)  \nonumber
\end{eqnarray}
    + Compute the local average estimate of the regression function at $x = 0.25$ and $x=0.75$ assuming that $h_{n} = 1/2$. 
    + Compute the k nearest neighbors estimate of the regression function at $x = 0.25$ and $x = 0.75$ assuming that $k = 2$.
    
\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

### The Nadaraya-Watson Estimator

* The Nadaraya-Watson estimator $\hat{m}_{h_{n}}^{NW}$ of the regression function with bandwidth $h_{n}$ is defined as
\begin{equation}
\hat{m}_{h_{n}}^{NW}(x) = \frac{ \sum_{i=1}^{n} Y_{i}K\Big( \frac{x - x_{i}}{ h_{n} }\Big) }{ \sum_{i=1}^{n} K\Big( \frac{x - x_{i}}{ h_{n} }\Big)  } \nonumber
\end{equation}

* The Nadaraya-Watson estimator has the same basic form as the local average estimator. We have just replaced the "box" function
$w(t)$ with the kernel function $K(t)$.

* You can think of $\hat{m}_{h_{n}}^{NW}(x)$ as a weighted average of the $Y_{i}$. 
That is, 
\begin{equation}
\hat{m}_{h_{n}}^{NW}(x) = \sum_{i=1}^{n} a_{i}(x) Y_{i} \nonumber
\end{equation}

* The bandwidth $h_{n}$ can also be referred to as the "smoothing parameter" since its value affects how smooth
the fitted regression curve appears.

* The weights $a_{1}(x), \ldots, a_{n}(x)$, in this case, are defined as 
\begin{equation}
a_{i}(x) = \frac{ K(\tfrac{x - x_{i}}{h_{n}})}{ \sum_{i=1}^{n} K(\tfrac{x - x_{i}}{ h_{n}}) } \nonumber
\end{equation}
So, we are using weights which are larger the closer you are to $x$. 

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* The Nadaraya-Watson estimator suffers from two main drawbacks. These are
**design bias** and **boundary bias**.

* Design bias refers to the effect of the spacing of the $x_{i}$ on the performance 
of the Nadaraya-Watson estimator.

* Boundary bias refers to the performance of the Nadaraya-Watson estimator
near the smallest and largest $x_{i}$.


\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* If we assume that the $x_{i}$ are random and have probability density $f_{X}(x)$, then it can be shown
that the mean-squared error of the Nadaraya-Watson estimator at a particular point $x$ has 
the following approximation
\begin{eqnarray}
\textrm{MSE}(x) &=& E\Big[ \{ m(x) - \hat{m}_{h_{n}}^{NW}(x) \}^{2} \Big] \nonumber \\
&\approx& \frac{h_{n}^{4}\mu_{2}^{2}(K)}{4}\Bigg\{ m''(x) + \frac{2m'(x)f_{X}'(x)}{f_{X}(x)} \Bigg\}^{2} + \frac{\sigma^{2}}{n h_{n} f_{X}(x) }, \nonumber
\end{eqnarray}
where $\mu_{2}(K) = \int_{-\infty}^{\infty} u^{2} K(u) du$ and $\kappa_{2}(K) = \int_{-\infty}^{\infty} K^{2}(u) du$. 

* The term $2m'(x)f_{X}'(x)/f_{X}(x)$ is referred to as the design bias. Notice that this should be zero if 
the $x_{i}$ are drawn from a Uniform distribution. In other words, if the $x_{i}$ are roughly equally spaced, then
the design bias should be small. 

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

**The Nadaraya-Watson estimator in R** 

* The Nadaraya-Watson estimator can be computed in `R` with the `ksmooth` function.

```r
ksmooth(x, y, kernel, bandwidth, x.points, ...)
```

*
    + `x` - vector of covariate values
    + `y` - vector of responses
    + `kernel` - choice of kernel function; default is `box`; use `normal` if you want a Gaussian kernel
    + `bandwidth` - value of the bandwidth; default is $0.05$
    + `x.points` - points at which to estimate the regression function; default is to use $n$ equally spaced points.

* The `x` vector from the fitted `ksmooth` object will be the vector of points at which the regression function
is estimated. The `y` vector from the fitted `ksmooth` object will be a vector containing the estimated values of the regression function.

* Note that the bandwidth used by this function for the Gaussian kernel is approximately $2.7$ times smaller than the bandwith in our definition of the Gaussian kernel.

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* If you wanted to write your own function that computed the Nadaraya-Watson estimate
at a vector of desired points $x.points = (t_{1}, \ldots, t_{q})$, you could use 
something like

```r
MyNWEst <- function(x, y, bandwidth, x.points) {
    q <- length(x.points) 
    nw.est <- numeric(q)
    for(k in 1:q) {
        ww <- dnorm(x.points[k], mean=x, sd=bandwidth) 
        nw.est[k] <- sum(ww*y)/sum(ww)
    }    
    return(nw.est)
}
```

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* To compute the Nadraya-Watson estimate at a set of equally spaced of points from $10$ to $25$
using bandwidth $0.5$ and plot the result, you could use the following code:

```r
tt <- seq(10, 25, by=.1)
bone.nwest <- ksmooth(x=bonedat$age, y=bonedat$spnbmd, kernel="normal", bandwidth=2.7*0.5, 
                      x.points=tt)

plot(bonedat$age, bonedat$spnbmd, las=1, ylab="Relative Change in Bone MD", 
     xlab="Age", main="Bone Data: Nadaraya-Watson Estimate with hn=0.5 and 
     Gaussian Kernel", type="n")
points(bonedat$age, bonedat$spnbmd, pch=16, cex=0.7)
lines(bone.nwest$x, bone.nwest$y, lwd=3, col="red")
```

![](11-kernel-regression-Latex_files/figure-latex/unnamed-chunk-9-1.pdf)<!-- --> 

```r
## Note that bone.nwest$x should equal tt
```



## Local Linear Regression

### Definition

* You can think of both the regressogram and the local average as methods
which fit local intercept models.

* For the regressogram, we fit an intercept model (that is a flat line curve) within 
a small bin that contains $x$. 

* For the local average estimator, we fit an intercept model for a small 
bin around $x$.

* The Nadaraya-Watson estimator can be thought of as just smoothing out the 
local intercept approach of the local average estimator. 


\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* Instead of fitting local intercept models, we could fit local
linear models that have an intercept and a slope term.

* To be specific, suppose we estimated the regression function at $x$
by fitting a linear model using only data from the $x_{i}$ that fell into the bin
$(x - h_{n}, x + h_{n})$.

* In this case, we would first fit the linear model 
$\hat{s}_{x}(x_{i}) = \hat{\beta}_{0x} + \hat{\beta}_{1x}(x_{i} - x)$ where $\hat{\beta}_{0x}$, $\hat{\beta}_{1x}$
solved the following local least-squares problem
\begin{equation}
\hat{\beta}_{0x}, \hat{\beta}_{1x} \textrm{ minimize: } \quad \sum_{i=1}^{n}\{ Y_{i} - \beta_{0x} - \beta_{1x}(x_{i} - x) \}^{2}I\big( x - h_{n} < x_{i} < x + h_{n} \big)  
(\#eq:simple-loclin-regression)
\end{equation}

* Then, we would estimate $m(x)$ by using the value of $\hat{s}_{x}$ at $x$. That is, $\hat{s}_{x}(x) = \hat{\beta}_{0x}$. 

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* **Local Linear Regression** uses the same idea as \@ref(eq:simple-loclin-regression), but replaces the indicator function $I( x - h_{n} < x_{i} < x + h_{n})$ with a smooth kernel function. 

* So, the local linear regression estimate of the regression function at $x$ is
\begin{eqnarray}
\hat{m}_{h_{n}}^{loclin}(x) &=& \hat{\beta}_{0x} \quad \textrm{ where } \nonumber \\
\hat{\beta}_{0x}, \hat{\beta}_{1x} &=&  \textrm{argmin}_{\beta_{0x},\beta_{1x}} \sum_{i=1}^{n}\{ Y_{i} - \beta_{0x} - \beta_{1x}(x_{i} - x) \}^{2}K\Big( \frac{x - x_{i}}{h_{n}} \Big) \nonumber
\end{eqnarray}

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* **Exercise 8.3** Suppose we define an estimator $\tilde{m}_{h_{n}}(x)$ of the regression function as
\begin{eqnarray}
\tilde{m}_{h_{n}}(x) &=& \hat{\beta}_{0x} \quad \textrm{ where } \nonumber \\
\hat{\beta}_{0x} &=&  \textrm{argmin}_{\beta_{0x}} \sum_{i=1}^{n}\{ Y_{i} - \beta_{0x} \}^{2}K\Big( \frac{x - x_{i}}{h_{n}} \Big) \nonumber
\end{eqnarray}
Show that $\tilde{m}_{h_{n}}(x) = \hat{m}_{h_{n}}^{NW}(x)$.

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

### Advantages of the Local Linear Estimator

* The local linear regression estimator can reduce the effects of design and boundary bias.

* If we write the local linear estimate at the point $x$ as
$\hat{m}_{h_{n}}^{loclin}(x) = \sum_{i=1}^{n} a_{i}^{h_{n}}(x)$, then the bias is appoximately 
\begin{equation}
E\{ \hat{m}_{h_{n}}^{loclin}(x) \} - m(x)
\approx m'(x)\sum_{i=1}^{n} (x_{i} - x)a_{i}^{h_{n}}(x) + \frac{m''(x)}{2}\sum_{i=1}^{n} (x_{i} - x)^{2}a_{i}^{h_{n}}(x) \nonumber 
\end{equation}


* For local linear regression, the term $m'(x)\sum_{i=1}^{n} (x_{i} - x)a_{i}^{h_{n}}(x)$ equals zero. If 
the weights $a_{i}^{h_{n}}(x)$ were the weights from the Nadaraya-Watson estimator, this term would not necessarily
equal zero.


* Also, the local linear estimator can help to reduce the boundary bias that arises
from asymmetry near the boundary (draw a picture).

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}


```r
MyLocLinear <- function(x, y, bandwidth, x.points) {
  q <- length(x.points) 
  loclin.est <- numeric(q)
  for(k in 1:q) {
    xtmp <- x - x.points[k]
    ww <- dnorm(xtmp, mean=0, sd=bandwidth) 
    loclin.est[k] <- lm(y ~ xtmp, weights=ww)$coef[1]
  }    
  return(loclin.est)
}
```

### Local Polynomial Regression


## Selecting the Bandwidth/Smoothing Parameter

### Representing in Linear Form

* Let $\mathbf{Y} = (Y_{1}, \ldots, Y_{n})$ and let $\hat{\mathbf{m}} = (\hat{m}(x_{1}), \ldots, \hat{m}(x_{n}))$ denote
the vector of "fitted values" from a vector of estimates of the regression function at $x_{1}, \ldots, x_{n}$.

* You can represent the fitted values all of the nonparametric estimators discussed thus far as
\begin{equation}
\hat{\mathbf{m}} = \mathbf{A}_{h_{n}}\mathbf{Y} \nonumber
\end{equation}
for an appropriately chosen $n \times n$ matrix $\mathbf{A}_{h_{n}}$.

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* For the local average estimator, we have $\hat{\mathbf{m}} = \mathbf{A}_{h_{n}}\mathbf{Y}$ where 
$\mathbf{A}_{h_{n}}$ is defined as
\begin{equation}
\mathbf{A}_{h_{n}}
= \begin{bmatrix} \frac{1}{n_{h_{n}}(x_{1})}I(x_{1} - h_{n} < x_{1} < x_{1} + h_{n}) & \ldots &  \frac{1}{n_{h_{n}}(x_{1})}I(x_{1} - h_{n} < x_{n} < x_{1} + h_{n}) \\
\frac{1}{n_{h_{n}}(x_{2})}I(x_{2} - h_{n} < x_{1} < x_{2} + h_{n}) & \ldots &  \frac{1}{n_{h_{n}}(x_{2})}I(x_{2} - h_{n} < x_{n} < x_{2} + h_{n}) \\
\vdots & \ddots & \vdots \\
\frac{1}{n_{h_{n}}(x_{n})}I(x_{n} - h_{n} < x_{1} < x_{n} + h_{n}) & \ldots &  \frac{1}{n_{h_{n}}(x_{n})}I(x_{n} - h_{n} < x_{n} < x_{n} + h_{n})
\end{bmatrix}
\end{equation}
where $n_{h_{n}}(x) = \sum_{i=1}^{n}I(x - h_{n} < x_{i} < x + h_{n})$.

* In other words, the $(i,j)$ element of $\mathbf{A}_{h_{n}}$ is $a_{i}(x_{j})$ where
\begin{equation}
a_{i}(x_{j}) = \frac{1}{n_{h_{n}}(x_{i})}I(x_{i} - h_{n} < x_{j} < x_{i} + h_{n})  \nonumber 
\end{equation}


\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* For the Nadaraya-Watson estimator, the $\mathbf{A}_{h_{n}}$ matrix is
\begin{equation}
\mathbf{A}_{h_{n}}
= \begin{bmatrix} \frac{1}{K_{h_{n}}(x_{1}, \cdot)  }K(0) &  \frac{1}{K_{h_{n}}(x_{1}, \cdot)  }K(\tfrac{x_{2} - x_{1}}{h_{n}}) & \ldots &  \frac{1}{K_{h_{n}}(x_{1}, \cdot)  }K(\tfrac{x_{n} - x_{1}}{h_{n}}) \\
\frac{1}{K_{h_{n}}(x_{2}, \cdot)  }K(\tfrac{x_{1} - x_{2}}{h_{n}}) &  \frac{1}{K_{h_{n}}(x_{2}, \cdot)  }K(0) & \ldots &  \frac{1}{K_{h_{n}}(x_{2}, \cdot)  }K(\tfrac{x_{n} - x_{2}}{h_{n}}) \\
\vdots & \vdots & \ddots & \vdots \\
\frac{1}{K_{h_{n}}(x_{n}, \cdot)  }K(\tfrac{x_{1}-x_{n}}{h_{n}}) &  \frac{1}{K_{h_{n}}(x_{n}, \cdot)  }K(\tfrac{x_{2} - x_{n}}{h_{n}}) & \ldots &  \frac{1}{K_{h_{n}}(x_{n}, \cdot)  }K(0)
\end{bmatrix} \nonumber
\end{equation}
where 
\begin{equation}
K_{h_{n}}(x_{i}, \cdot) = \sum_{j=1}^{n}K\Big( \frac{x_{i} - x_{j}}{h_{n}}  \Big) \nonumber
\end{equation}

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* For the local linear regression estimator, 

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* In "classic" linear regression where you would try to fit the straight line model 
$Y_{i} = \beta_{0} + \beta_{1}x_{i} + \varepsilon_{i}$, the vector of fitted values would be
\begin{equation}
\hat{\mathbf{m}} = (\hat{m}(x_{1}), \ldots, \hat{m}(x_{n}))
= (\hat{\beta}_{0} + \hat{\beta}_{1}x_{1}, \ldots, \hat{\beta}_{0} + \hat{\beta}_{1}x_{n})  \nonumber
\end{equation}

* In this case, you can represent $\hat{\mathbf{m}}$ as
\begin{equation}
\hat{\mathbf{m}} = \mathbf{X}(\mathbf{X}^{T}\mathbf{X})^{-1}\mathbf{X}^{T}\mathbf{Y}  \nonumber
\end{equation}
where $\mathbf{X}$ is the following $n \times 2$ "design" matrix
\begin{equation}
\mathbf{X} = \begin{bmatrix} 1 & x_{1} \\ 1 & x_{2} \\ \vdots & \vdots \\ 1 & x_{n} \end{bmatrix} \nonumber
\end{equation}


### The Cp Statistic

* **Theorem:** If a random vector $\mathbf{Z}$ has mean vector $\mathbf{\mu}$ and covariance matrix $\mathbf{\Sigma}$,
then
\begin{equation}
E\{ \mathbf{Z}^{T}\mathbf{Z} \} = E\Big\{ \sum_{i=1}^{n} Z_{i}^{2} \Big\} = \mathbf{\mu}^{T}\mathbf{\mu} + \textrm{tr}( \mathbf{\Sigma} ) \nonumber
\end{equation}

* Notice that the vector $\mathbf{m} - \mathbf{A}_{h_{n}}\mathbf{Y}$ has 
\begin{equation}
E( \mathbf{m} - \mathbf{A}_{h_{n}}\mathbf{Y} ) = (\mathbf{I} - \mathbf{A}_{h_{n}})\mathbf{m}  \qquad \qquad \textrm{Var}(\mathbf{m} - \mathbf{A}_{h_{n}}\mathbf{Y}) = \sigma^{2}\mathbf{A}_{h_{n}}\mathbf{A}_{h_{n}}^{T}
\end{equation}

* Also, the vector $\mathbf{Y} - \mathbf{A}_{h_{n}}\mathbf{Y} = (\mathbf{I} - \mathbf{A}_{h_{n}})\mathbf{Y}$ has
\begin{equation}
E\{ (\mathbf{I} - \mathbf{A}_{h_{n}})\mathbf{Y} \} = (\mathbf{I} - \mathbf{A}_{h_{n}})\mathbf{m}  \qquad \qquad \textrm{Var}\{ (\mathbf{I} - \mathbf{A}_{h_{n}})\mathbf{Y}) = \sigma^{2} (\mathbf{I} - \mathbf{A}_{h_{n}})(\mathbf{I} - \mathbf{A}_{h_{n}})^{T} \nonumber
\end{equation}

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* Ideally, we would like to choose the smoothing parameter $h_{n}$ to minimize the following mean averaged squared error
\begin{eqnarray}
\textrm{MASE}(h_{n}) &=& \frac{1}{n}\sum_{i=1}^{n} E\Big[ \{ m(x_{i}) - \hat{m}(x_{i}) \}^{2} \Big]  \nonumber \\
&=& E\{ \frac{1}{n}[\mathbf{m} - \mathbf{A}_{h_{n}}\mathbf{Y}]^{T}[\mathbf{m} - \mathbf{A}_{h_{n}}\mathbf{Y}] \} \nonumber
\end{eqnarray}

* If we apply the above Theorem to the vector $\mathbf{m} - \mathbf{A}_{h_{n}}\mathbf{Y}$, we can notice that
\begin{eqnarray}
\textrm{MASE}( h_{n} ) &=& E\{ \frac{1}{n}( \mathbf{m} - \mathbf{A}_{h_{n}}\mathbf{Y} )^{T}(\mathbf{m} - \mathbf{A}_{h_{n}}\mathbf{Y}) \}  \nonumber \\
&=& \frac{1}{n}[(\mathbf{I} - \mathbf{A}_{h_{n}})\mathbf{m}]^{T}[(\mathbf{I} - \mathbf{A}_{h_{n}})\mathbf{m}] + \frac{\sigma^{2}}{n}\textrm{tr}(\mathbf{A}_{h_{n}}\mathbf{A}_{h_{n}}^{T})
\end{eqnarray}

* Now, using the mean and covariance matrix for $(\mathbf{I} - \mathbf{A}_{h_{n}})\mathbf{Y}$, we also have that
\begin{eqnarray}
&& E\{\frac{1}{n}[(\mathbf{I} - \mathbf{A}_{h_{n}})\mathbf{Y}]^{T}[(\mathbf{I} - \mathbf{A}_{h_{n}})\mathbf{Y}] \} \nonumber \\
&=& \frac{1}{n}[(\mathbf{I} - \mathbf{A}_{h_{n}})\mathbf{m}]^{T}[(\mathbf{I} - \mathbf{A}_{h_{n}})\mathbf{m}] + \frac{\sigma^{2}}{n}\textrm{tr}\{ (\mathbf{I} - \mathbf{A}_{h_{n}})(\mathbf{I} - \mathbf{A}_{h_{n}})^{T} \} \nonumber \\
&=& \textrm{MASE}(h_{n}) + \sigma^{2} - \frac{2\sigma^{2}}{n}\textrm{tr}( \mathbf{A}_{h_{n}})
\end{eqnarray}

* So, if $\sigma^{2}$ is known, then $\widehat{\textrm{MASE}}( h_{n} )$ is an unbiased estimate of $\textrm{MASE}( h_{n} )$
\begin{eqnarray}
\widehat{\textrm{MASE}}( h_{n} ) &=& \frac{1}{n}[(\mathbf{I} - \mathbf{A}_{h_{n}})\mathbf{Y}]^{T}[(\mathbf{I} - \mathbf{A}_{h_{n}})\mathbf{Y}] - \sigma^{2} + \frac{2\sigma^{2}}{n}\textrm{tr}( \mathbf{A}_{h_{n}}) \nonumber \\
&=& \frac{1}{n}\sum_{i=1}^{n} \{ Y_{i} - \hat{m}_{h_{n}}(x_{i}) \}^{2} - \sigma^{2} + \frac{2\sigma^{2}}{n}\textrm{tr}( \mathbf{A}_{h_{n}})  \nonumber
\end{eqnarray}

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* The predictive mean averaged squared error (PAMSE) is defined as
\begin{equation}
PAMSE(h_{n}) = E\Big[ \frac{1}{n} \sum_{i=1}^{n} \{ Y_{i}' -  \hat{m}_{h_{n}}(x_{i}) \}^{2} \Big]
\end{equation}
where $Y_{i}'$ is a "future" independent observation that has the same covariate as $Y_{i}$. 

* We assume that $Y_{i}' = m(x_{i}) + \varepsilon_{i}'$ where $\varepsilon_{i}'$ is independent of $\varepsilon_{i}$.

* So, 
\begin{eqnarray}
\textrm{PAMSE}(h_{n}) &=& E\Big[ \frac{1}{n} \sum_{i=1}^{n} \{ m(x_{i}) -  \hat{m}_{h_{n}}(x_{i}) + \varepsilon_{i}' \}^{2} \Big] 
\nonumber \\
&=& E\Big[ \frac{1}{n} \sum_{i=1}^{n} \{ m(x_{i}) -  \hat{m}_{h_{n}}(x_{i}) \}^{2} \Big]  +  E\Big[ \frac{1}{n} \sum_{i=1}^{n} (\varepsilon_{i}')^{2} \Big]  \nonumber \\
&=& \textrm{AMSE}( h_{n} ) + \sigma^{2} \nonumber
\end{eqnarray}


\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* The $C_{p}$ statistic is based on the idea that, if $\sigma^{2}$ was known, then the following quantity
would be an unbiased estimate of $\textrm{PAMSE}( h_{n} )$:
\begin{equation}
\frac{1}{n}\sum_{i=1}^{n} \{ Y_{i} - \hat{m}_{h_{n}}(x_{i}) \}^{2}  + \frac{2\sigma^{2}}{n}\textrm{tr}( \mathbf{A}_{h_{n}}) \nonumber
\end{equation}

* The $C_{p}$ statistic formula is obtained by plugging in an estimate $\hat{\sigma}^{2}$ of the residual variance into the above
formula:
\begin{equation}
C_{p} = \frac{1}{n}\sum_{i=1}^{n} \{ Y_{i} - \hat{m}_{h_{n}}(x_{i}) \}^{2}  + \frac{2\hat{\sigma}^{2}}{n}\textrm{tr}( \mathbf{A}_{h_{n}}) \nonumber
\end{equation}



### Leave-one-out Cross Validation

* Similar to the way we defined a leave-on-out density estimate in Chapter 8, we will define
the leave-one-out estimate of the regression function at $x$ as:
\begin{equation}
\hat{m}_{h_{n},-i}(x) - \textrm{ estimate of $m(x)$ found by using all data except $(Y_{i}, x_{i})$.}
\end{equation}


\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* Because we are assuming that $\hat{m}_{h_{n}}(x)$ can be represented as a linear combination of the responses
\begin{equation}
\hat{m}_{h_{n}}(x) = \sum_{j=1}^{n} a_{j}^{h_{n}}(x)Y_{j}
\end{equation}
the leave-one-out estimate can be expressed as
\begin{equation}
\hat{m}_{h_{n}, -i}(x) = \sum_{j=1}^{n} a_{j,-i}^{h_{n}}(x)Y_{j}
\end{equation}
where 
\begin{equation}
a_{j,-i}^{h_{n}}(x) = 
\end{equation}

\begin{center}
\rule{\textwidth}{.05cm}
\end{center}

* The predictive risk estimate is



## Additional Reading

* Additional reading which covers the material discussed in this chapter includes:
    + Chapter 4 from @hardle2012
    + Chapter 5 from @wasserman2006
