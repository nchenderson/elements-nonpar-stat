# (PART) Nonparametric Estimation {-} 

# The Empirical Distribution Function {#edf}


## Definition and Basic Properties

* Every random variable has a cumulative distribution function (cdf).

* The cdf of a random variable $X$ is defined as
\begin{equation}
F(t) = P( X \leq t)
\end{equation}

* The empirical distribution function or empirical cumulative distribution function (ecdf)
estimates $F(t)$ by computing the proportion of observations which are less than or equal 
to $t$. 

* For i.i.d. random variables $X_{1}, \ldots, X_{n}$ with cdf $F$, the empirical distribution function
is defined as
\begin{equation}
\hat{F}_{n}(t) = \frac{1}{n}\sum_{i=1}^{n} I( X_{i} \leq t) \nonumber
\end{equation}

* Note that the empirical distribution function can be computed for any type
of data without making any assumptions about the distribution from which 
the data arose.

* The only assumption we are making is that $X_{1}, \ldots, X_{n}$
constitute an i.i.d. sample from some common distribution function $F$.

---


## Confidence intervals for F(t)

* For a fixed value of $t$, the distribution of $\hat{F}_{n}(t)$ is
\begin{equation}
n \hat{F}_{n}(t) \sim \textrm{Binomial}\big( n, F(t) \big)
\end{equation}
(why?)

* For a fixed $t$, $n\hat{F}_{n}(t)$ is the sum of $n$ independent
Bernoulli random variables $W_{1}^{t}, \ldots, W_{n}^{t}$
\begin{equation}
\hat{F}_{n}(t) = \sum_{i=1}^{n} W_{i}^{t} = \sum_{i=1}^{n} I( X_{i} \leq t)
\end{equation}

* The probability that $W_{i}^{t} = 1$ is
\begin{equation}
P( W_{i}^{t} = 1) = P(X_{i} \leq t) = F(t) \nonumber
\end{equation}

---



## The Empirical Distribution Function in R

* We will see how to work with empirical distribution functions in **R** by using the "kidney data."


```r
kidney <- read.table("https://web.stanford.edu/~hastie/CASI_files/DATA/kidney.txt", 
                     header=TRUE)
head(kidney)
```

```
##   age   tot
## 1  18  2.44
## 2  19  3.86
## 3  19 -1.22
## 4  20  2.30
## 5  21  0.98
## 6  21 -0.50
```

* The **ecdf** function is the main function which computes the empirical distribution function
in **R**

* The **ecdf** function will create an **ecdf** object. To create an ecdf object
for the kidney totals, use the following code:

```r
kidney.Fhat <- ecdf(kidney$tot)
```

* You can plot the ecdf for the kidney totals by just calling **plot(ecdf)**

```r
plot(kidney.Fhat)
```

<img src="07-empiricalcdf_files/figure-html/unnamed-chunk-3-1.png" width="672" />

* You can plot the ecdf for the kidney totals by just calling **plot(ecdf)**

```r
plot(kidney.Fhat, verticals=TRUE)
```

<img src="07-empiricalcdf_files/figure-html/unnamed-chunk-4-1.png" width="672" />

## The Kolmogorov-Smirnov Test






## The empirical distribution function and statistical functionals

