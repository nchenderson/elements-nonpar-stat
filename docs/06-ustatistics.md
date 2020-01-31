
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
U &=& \frac{1}{{n \choose 2}} \sum_{p \in \mathcal{P}_{n,2}} h(X_{p_{1}}, X_{p_{2}})
= \frac{2}{n(n-1)}\sum_{i=1}^{n} \sum_{j=i+1}^{n} \frac{1}{2} (X_{i} - X_{j})^{2}  \nonumber \\
&=& \frac{2}{n(n-1)}\sum_{i=1}^{n} \sum_{j=1}^{n} (X_{i} - X_{j})^{2} \nonumber \\
&=& \frac{2}{n(n-1)}\sum_{i=1}^{n} \sum_{j=1}^{n} \{ (X_{i} - \bar{X})^{2} - 2(X_{i} - \bar{X})(X_{j} - \bar{X}) + (X_{j} - \bar{X})^{2} \} \nonumber \\
&=& \frac{2}{n(n-1)}\sum_{i=1}^{n} n(X_{i} - \bar{X})^{2} + \frac{2}{n(n-1)}\sum_{j=1}^{n} n(X_{j} - \bar{X})^{2} \} \nonumber \\
&=& \frac{1}{n-1}\sum_{i=1}^{n} (X_{i} - \bar{X})^{2}
\end{eqnarray}

### Example 3: Gini's Mean Difference

* Gini's mean difference statistic is defined as
\begin{equation}
\frac{1}{{n \choose 2}} \sum_{i=1}^{n}\sum_{j=i+1} | X_{i} - X_{j} | \nonumber
\end{equation}

* This is a U-statistic with kernel
\begin{equation}
h(X_{1},X_{2}) = | X_{1} - X_{2} |
\end{equation}

* The parameter that we are estimating with Gini's mean difference statistic is:
\begin{equation}
\theta = E\Big\{ \Big| X_{1} - X_{2} \Big|  \Big\}
\end{equation}

* We can interpret $\theta$ in the following way: If we draw 
two observations at random from our population, $\theta$ 
represents the expected absolute difference between these
two observations.

### Example 4: Wilcoxon Signed Rank Statistic


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
U = \frac{1}{mn}\sum_{i=1}^{n}\sum_{j=1}^{m} I( X_{i} \leq Y_{j})
\end{equation}

* This is a U-statistic with kernel $h(x, y) = I(x \leq y)$.

* The statistic $mn U$ is known as the **Mann-Whitney** statistic.

* Relation to Wilcoxon rank sum statistic?


## Measures of Association

### Spearman's Rank Correlation


### Kendall's tau


### Distance Covariance


