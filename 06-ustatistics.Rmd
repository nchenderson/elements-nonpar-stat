
# U-Statistics {#ustat}

## Definition

* Suppose we have observations $X_{1}, \ldots, X_{n}$.

* U-statistics are a family of statistics used to estimate quantities 
that can be written as
\begin{equation}
\theta = E\Big\{ h(X_{1}, \ldots, X_{r})  \Big\}
(\#eq:ustat-parameter)
\end{equation}

* The U-statistic $U$ which estimates \@ref(eq:ustat-parameter) is given by the following formula:
\begin{equation}
U = \frac{1}{{n \choose r}} \sum_{p \in \mathcal{P}_{n,r}} h(X_{p_{1}}, \ldots, X_{p_{r}})
(\#eq:ustat-definition)
\end{equation}

* The function $h$ is usually called the **kernel** of the U-statistic.

* The integer $r$ is called the **order** of the U-statistic. Typically,
$r = 2$, or $r = 3$ at most.



## Examples
 
* A wide range of well-known statistics can be represented as U-statistics.


### Example 1: The Sample Mean

* The sample mean is actually an example of a U-statistic with $r = 1$.

* Choosing $h(x) = x$ means that the corresponding U-statistic is
\begin{equation}
U = \frac{1}{n} \sum_{i=1}^{n} X_{i}
\end{equation}

### Example 2: The Sample Variance

* The sample variance is actually another example of a U-statistic. 
In this case, $r = 2$.

* To show why this is the case, choose the kernel $h(x_{1}, x_{2})$ to be
\begin{equation}
h(x_{1}, x_{2}) = \frac{1}{2}(x_{1} - x_{2})^{2}
\end{equation}

* The expectation of this kernel is $\sigma^{2} = E\{ h(X_{1}, X_{2}) \}$ because
\begin{eqnarray}
E\{ h(X_{1}, X_{2}) \} &=& \frac{1}{2}\Big[ E(X_{1}^{2}) - 2E(X_{1})E(X_{2})  + E(X_{2}^{2}) \Big] \nonumber \\
&=& \frac{1}{2}\Big[ \sigma^{2} + \mu^{2} - 2\mu^{2}  + \sigma^{2} + \mu^{2} \Big] \nonumber \\
&=& \sigma^{2}
\end{eqnarray}

* Also, by using formula \@ref(eq:ustat-definition), this choice of kernel generates the sample variance at its U-statistic:
\begin{eqnarray}
U_{var} &=& \frac{1}{{n \choose 2}} \sum_{p \in \mathcal{P}_{n,2}} h(X_{p_{1}}, X_{p_{2}})
= \frac{2}{n(n-1)}\sum_{i=1}^{n} \sum_{j=i+1}^{n} \frac{1}{2} (X_{i} - X_{j})^{2}  \nonumber \\
&=& \frac{2}{n(n-1)}\sum_{i=1}^{n} \sum_{j=1}^{n} (X_{i} - X_{j})^{2} \nonumber \\
&=& \frac{2}{n(n-1)}\sum_{i=1}^{n} \sum_{j=1}^{n} \{ (X_{i} - \bar{X})^{2} - 2(X_{i} - \bar{X})(X_{j} - \bar{X}) + (X_{j} - \bar{X})^{2} \} \nonumber \\
&=& \frac{2}{n(n-1)}\sum_{i=1}^{n} n(X_{i} - \bar{X})^{2} + \frac{2}{n(n-1)}\sum_{j=1}^{n} n(X_{j} - \bar{X})^{2} \} \nonumber \\
&=& \frac{1}{n-1}\sum_{i=1}^{n} (X_{i} - \bar{X})^{2}
\end{eqnarray}

### Example 3: Gini's Mean Difference

* Gini's mean difference statistic is defined as
\begin{equation}
U_{G} = \frac{1}{{n \choose 2}} \sum_{i=1}^{n}\sum_{j=i+1} | X_{i} - X_{j} | \nonumber
\end{equation}

* This is a U-statistic with kernel
\begin{equation}
h(X_{1},X_{2}) = | X_{1} - X_{2} |
\end{equation}

* The parameter that we are estimating with Gini's mean difference statistic is:
\begin{equation}
\theta_{G} = E\Big\{ \Big| X_{1} - X_{2} \Big|  \Big\}
\end{equation}

* Gini's mean difference parameter $\theta_{G}$ can be interpreted in the following way: If we draw 
two observations at random from our population, $\theta_{G}$ 
represents the expected absolute difference between these
two observations.

* The Gini coefficient $\theta_{Gc}$ is a popular measure of inequality. It is related
to the mean difference parameter via
\begin{equation}
\theta_{Gc} = \frac{ \theta_{G}}{ 2\mu },
\end{equation}
where $\mu = E( X_{i} )$.

---

* **Exercise 6.1**. Compute the Gini coefficient $\theta_{Gc}$ when it is assumed that
    + $X_{i} \sim \textrm{Normal}( \mu, \sigma^{2})$, for $\mu > 0$. 
    + $X_{i} \sim \textrm{Exponential}(\lambda)$, (**Hint**: The difference between two independent Exponential random variables has a Laplace distribution).

---

### Example 4: Wilcoxon Signed Rank Statistic

* The Wilcoxon signed rank test statistic has some relation to the following U statistic
\begin{equation}
U_{WS} = \frac{2}{n(n-1)}\sum_{i=1}^{n}\sum_{j=i+1}^{n} I\Big( X_{i} + X_{j} > 0 \Big)
\end{equation}

* Hence, the Wilcoxon signed ran test statistic can be interpreted as an estimate of
the following parameter
\begin{equation}
\theta_{WS} = P\Big( X_{i} + X_{j} > 0  \Big) = P\Big( X_{i} > -X_{j} \Big)
\end{equation}

* If the distribution of $X_{i}$ is symmetric around $0$, $\theta_{WS}$ will be equal 
to $1/2$.

* Recall that the Wilcoxon signed rank test is designed to detect
distributions which are not symmetric around $0$.

* $U_{WS} = a T_{n}$, where $T_{n}$ is the signed rank statistic that we defined in 
Section 3.3.

## U-statistics for Two-Sample Problems

* In two-sample problems, we have data from two groups which
we label $X_{1}, \ldots, X_{n}$ and $Y_{1}, \ldots, Y_{m}$

* A U-statistic with order $(r,s)$ for a two-sample problem is
\begin{equation}
U = \frac{1}{{n \choose r}}\frac{1}{{m \choose s}} \sum \sum h(X_{p_{1}}, \ldots, X_{p_{r}}, Y_{q_{1}}, \ldots, Y_{q_{s}})
\end{equation}


### The Mann-Whitney Statistic

* Consider the following U-statistic 
\begin{equation}
U_{MW} = \frac{1}{mn}\sum_{i=1}^{n}\sum_{j=1}^{m} I( X_{i} \geq Y_{j})
\end{equation}

* This is a U-statistic of order $(1,1)$ with kernel $h(x, y) = I(x \geq y)$.

* The U-statistic $U_{MW}$ can be thought of as an estimate of the following
parameter
\begin{equation}
\theta_{MW} = P\Big( X_{i} \geq Y_{j} \Big)
\end{equation}

* If both $X_{i}$ and $Y_{j}$ have the same distribution, then
$\theta_{MW}$ should equal $1/2$.

---


* The statistic $mn U_{MW}$ is known as the **Mann-Whitney** statistic.

* The Mann-Whitney statistic has a close relation to the Wilcoxon 
rank sum statistic $W$ that we defined in Section 3.2:
\begin{eqnarray}
mn U_{MW} &=& 
\sum_{i=1}^{n}\sum_{j=1}^{m} I( X_{i} \geq Y_{j}) \nonumber \\
&=& \sum_{i=1}^{n}\Big[ \sum_{j=1}^{m} I( X_{i} \geq Y_{j}) +
\sum_{k=1}^{n} I( X_{i} \geq X_{k}) \Big] -
\sum_{i=1}^{n}\sum_{k=1}^{n} I( X_{i} \geq X_{k}) \nonumber \\
&=& \sum_{i=1}^{n} R_{i}(\mathbf{Z}) -
\sum_{i=1}^{n} R_{i}( \mathbf{X} ) \nonumber \\
&=& W - \frac{n(n+1)}{2} \nonumber
\end{eqnarray}

* In other words, the Mann-Whitney statistic is equal to 
the WRS statistic minus a constant term.

* In the above derivation, we are using $\mathbf{Z}$ as
the pooled-data vector $\mathbf{Z} = (X_{1}, \ldots, X_{n}, Y_{1}, \ldots, Y_{m})$.

* Also, the above derivation assumes no ties so that 
$\sum_{i=1}^{n} R_{i}( \mathbf{X} ) = n(n+1)/2$.

---

* Comment on implications for alternative in the WRS test here.

## Measures of Association

* Many important measures of association are also examples of U-statistics.

* For measures of association, we have observations on $n$ pairs of variables
\begin{equation}
(X_{1}, Y_{1}), \ldots, (X_{n}, Y_{n}), \nonumber 
\end{equation}
and our goal is to report some measure which quantifies the relationship between these
two variables.

* In this context, we will think about U-statistics having the form
\begin{equation}
U = \frac{1}{{n \choose 2}} \sum_{p \in \mathcal{P}_{n,2} } h\Bigg( \begin{bmatrix} X_{p_{1}} \\ Y_{p_{1}} \end{bmatrix},
\begin{bmatrix} X_{p_{2}} \\ Y_{p_{2}} \end{bmatrix} \Bigg)
\end{equation}


### Spearman's Rank Correlation


### Kendall's tau

* Kendall's $\tau$-statistic $U_{\tau}$ is given by
\begin{eqnarray}
U_{\tau} &=& \frac{2}{n(n-1)}\sum_{i=1}^{n}\sum_{j=1}^{n} I\Big\{ (X_{j} - X_{i})(Y_{j} - Y_{i}) > 0 \Big\}  - 1 \nonumber \\
&=& \frac{4}{n(n-1)}\sum_{i=1}^{n}\sum_{j=i+1}^{n} I\Big\{ (X_{j} - X_{i})(Y_{j} - Y_{i}) > 0 \Big\}  - 1
\end{eqnarray}

* Note that $U_{\tau}$ is a U-statistic of order $2$ with kernel
\begin{equation}
h\Bigg( \begin{bmatrix} X_{1} \\ Y_{1} \end{bmatrix},
\begin{bmatrix} X_{2} \\ Y_{2} \end{bmatrix} \Bigg)
= 2 \times I\Big\{ (X_{2} - X_{1})(Y_{2} - Y_{1}) > 0 \Big\}  - 
I\Big\{ X_{2} \neq X_{1} \Big\}I\Big\{ Y_{2} \neq Y_{1} \Big\}
\end{equation}

* Assuming the probability of ties is zero, Kendall's $\tau$ can be thought of as an estimate of the following quantity
\begin{equation}
2 P\Big\{ (X_{j} - X_{i})(Y_{j} - Y_{i}) > 0 \Big\} - 1
\end{equation}

* If $X_{i}$ and $Y_{i}$ are idependent, Kendall's $\tau$ will be equal to zero (why?).

* Also, Kendall's $\tau$ must be in between $-1$ and $1$.

---

* In the context of computing $U_{\tau}$, pairs of observations $(X_{i}, Y_{i})$ and $(X_{j}, Y_{j})$ are said to be **concordant** if 
the sign of $X_{j} - X_{i}$ agrees with the sign of $Y_{j} - Y_{i}$.

* If the sign of $X_{j} - X_{i}$ and $Y_{j} - Y_{i}$ do **not** agree, then the
pairs $(X_{i}, Y_{i})$ and $(X_{j}, Y_{j})$ are said to be **discordant**.

* If either $X_{j}=X_{i}$ or $Y_{j}=Y_{i}$, then the pairs $(X_{i}, Y_{i})$ and $(X_{j}, Y_{j})$
are neither concordant or discordant.

* Let us defined
\begin{eqnarray}
n_{c} &=& \textrm{ the number of concordant pairs} \nonumber \\
n_{d} &=& \textrm{ the number of discordant pairs} \nonumber \\
n_{n} &=& \textrm{ number of pairs which are neither}  \nonumber
\end{eqnarray}

* Then, 
\begin{equation}
n_{c} + n_{d} + n_{n} = {n \choose 2} = \frac{n(n-1)}{2}
\end{equation}

* If we assume that there are no ties (i.e., $n_{n} = 0$), then $U_{\tau}$
can be written as
\begin{equation}
U_{\tau} = \frac{4n_{c}}{n(n-1)} - 1
= \frac{2n_{c} + 2n_{c} - n(n - 1)}{n(n-1)}
= \frac{2n_{c} - 2n_{d} }{n(n-1)} 
= \frac{2(n_{c} - n_{d})}{n(n-1)} \nonumber
\end{equation}

* Under independence, the number of concordant and discordant pairs should be roughly equal.

---


### Distance Covariance

* Define $a_{ij}$ and $b_{ij}$ as
\begin{eqnarray}
a_{ij} &=& | X_{i} - X_{j}|  \nonumber \\
b_{ij} &=& | Y_{i} - Y_{j}|  \nonumber
\end{eqnarray}



