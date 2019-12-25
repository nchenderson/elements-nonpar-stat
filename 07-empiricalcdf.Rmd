# (PART) Nonparametric Estimation {-} 

# The Empirical Distribution Function {#regression}

## Empirical Distribution Functions

### Definition and Basic Properties

* Every random variable has a cumulative distribution function (cdf).
* The cdf of a random variable $X$ is defined as
\begin{equation}
F(t) = P( X \leq t)
\end{equation}
* 

* The empirical distribution function or empirical cumulative distribution function (ecdf)
estimates $F(t)$ by just finding the proportion of observations which are less than or equal 
to $t$. 
* For i.i.d. random variables $X_{1}, \ldots, X_{n}$, the empirical distribution function
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


### Confidence intervals for Fhat

## The Empirical Distribution Function in R







## The empirical distribution functino and statistical functionals

