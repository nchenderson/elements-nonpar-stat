
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

* Also, this choice of kernel generates the sample variance because
\begin{equation}
U = \frac{1}{{n \choose 2}} \sum_{p \in \mathcal{P}_{n,2}} h(X_{p_{1}}, X_{p_{2}})
= \frac{1}{2}\sum_{i=1}^{n} \sum_{j=i+1}^{n} (X_{i} - X_{j})^{2}
\end{equation}

### Example 3: Gini's Mean Difference

* The parameter of interest is
\begin{equation}
\theta = E\Big\{ \Big| X_{1} - X_{2} \Big|  \Big\}
\end{equation}

## U-statistics for One and Two-Sample Problems

### The Mann-Whitney Statistic

## Measures of Association

### Spearman's Rank Correlation


### Kendall's tau


### Distance Covariance


