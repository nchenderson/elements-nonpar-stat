# (PART) Nonparametric Testing {-} 

# Rank and Sign Statistics {#rank-tests}

## Introduction

Start with t-test example, difference in means is sufficient for superiority

Give example of type of tests we are interested in.

Why ranks and why nonparametric testing?

(reduce influence of outliers)

## Ranks

### Definition

* Suppose we have $n$ observations $X_{1}, \ldots, X_{n}$. The **rank** of the $i^{th}$ observation $R_{i}$ is defined as
\begin{equation}
R_{i} = R(X_{i}) = \sum_{j=1}^{n} I( X_{j} \leq X_{i}) 
(\#eq:rankdef)
\end{equation}
where
\begin{equation}
I(X_{j} \leq X_{i}) 
= \begin{cases}
1 & \text{ if } X_{j} \leq X_{i} \\
0 & \text{ if } X_{j} > X_{i}
\end{cases}
\end{equation}
* The largest observation has a rank of $n$.
* The smallest observation has a rank of $1$ (if there are no ties).


```r
x <- c(3, 7, 1, 12, 6)  ## 5 observations
rank(x)
```

```
## [1] 2 4 1 5 3
```

* In the definition of ranks shown in \@ref(eq:rankdef), tied observations
receive their maximum possible rank. 
* For example, suppose that $(X_{1}, X_{2}, X_{3}, X_{4}) = (0, 1, 1, 2)$. 
In this case, one could argue whether both observations 2 and 3 should be ranked
$2^{nd}$ or $3^{rd}$ while observations $1$ and $4$ should unambiguously receive
ranks of $1$ and $4$ respectively.
* Under definition \@rank(eq:rankdef), both observations $2$ and $3$ receive a rank of $3$.

* In **R**, such handling of ties is done using the **ties.method = "max"** argument

```r
x <- c(0, 1, 1, 2)  
rank(x, ties.method="max")
```

```
## [1] 1 3 3 4
```
* The default in **R** is to replace the ranks of tied observations with their "average" rank

```r
x <- c(0, 1, 1, 2)  
rank(x)
```

```
## [1] 1.0 2.5 2.5 4.0
```

### Properties of Ranks
Suppose $(X_{1}, \ldots, X_{n})$ is random sample from a continuous distribution $F$ (so that the probability
of ties is zero). Then, the following properties hold for the associated ranks $R_{1}, \ldots, R_{n}$.

* Each $R_{i}$ follows a discrete uniform distribution
\begin{equation}
P(R_{i} = j) = 1/n, \quad \text{for any } j = 1, \ldots,n.
\end{equation}
* The expectation of $R_{i}$ is
\begin{equation}
E( R_{i} ) = \sum_{j=1}^{n} j P(R_{i} = j) = \frac{1}{n}\sum_{j=1}^{n} j = \frac{(n+1)}{2}
\end{equation}
* The variance of $R_{i}$ is
\begin{equation}
\text{Var}( R_{i} ) = E( R_{i}^{2} ) - E(R_{i})^{2}
= \frac{1}{n}\sum_{j=1}^{n} j^{2}  - \Big( \frac{n+1}{2} \Big)^{2}
= \frac{ n^{2} - 1}{12}
\end{equation}
* The random variables $R_{1}, \ldots, R_{n}$ are **not** independent (why?). However,
the vector $\mathbf{R}_{n} = (R_{1}, \ldots, R_{n})$ is uniformly distributed
on the set of $n!$ permutations of $(1,2,\ldots,n)$.


## Two-Sample Tests

### The Wilcoxon Rank Sum (WRS) Test

#### Goal of the Test

* The Wilcoxon Rank Sum (WRS) test (sometimes referred to as the Wilcoxon-Mann-Whitney test) is a popular,
rank-based two-sample test.

* The WRS test is used to test whether or not observations from one group tend to be larger (or smaller) than observations
from the other group. 

* Suppose we have observations from two groups: $X_{1}, \ldots, X_{n} \sim F_{X}$ and $Y_{1}, \ldots, Y_{m} \sim F_{Y}$.

* Roughly speaking, the WRS tests the following hypothesis
\begin{eqnarray}
H_{0}: & & F_{X} = F_{Y} \quad \textrm{ versus } \nonumber \\
H_{A}: & & \textrm{Observations from } F_{X} \textrm{ tend to be larger than observations from } F_{Y} \nonumber 
\end{eqnarray}

---

* What is meant by "tend to be larger" in the alternative hypothesis?

* Two common ways of stating the alternative hypothesis for the WRS include
    1. The stochastic dominance alternative
\begin{eqnarray}
H_{0}: & & F_{X} = F_{Y} \quad \textrm{ versus } \nonumber \\
H_{A}: & & F_{X} \textrm{ is stochastically larger than } F_{Y} \nonumber 
\end{eqnarray}
    2. The "shift" alternative
\begin{eqnarray}
H_{0}: & & F_{X} = F_{Y} \quad \textrm{ versus } \nonumber \\
H_{A}: & & F_{X}(t) = F_{Y}(t - \Delta), \Delta > 0.
\end{eqnarray}
* A distribution function $F_{X}$ is said to be stochastically larger than
$F_{Y}$ if $F_{X}(t) \geq F_{Y}(t)$ for all $t$ with $F_{X}(t) > F_{Y}(t)$
for at least one value of $t$.

* Note that the "shift alternative" implies stochastic dominance.

---

* Why do we need to specify an alternative?

#### Definition of the WRS Test Statistic

* The WRS test statistic is based on computing the sum of ranks (ranks based on the pooled sample)
in one group.

* If observations from group 1 tend to be larger than those from group 2, the average rank from group 1 should exceed the
average rank from group 2. 

* A sufficiently large value of the average rank from group 1 will allow to reject $H_{0}$ 
in favor of $H_{A}$.

---

* Let us define the pooled sample data vector $\mathbf{Z}$ as 
\begin{equation}
\mathbf{Z} = (X_{1}, \ldots, X_{n}, Y_{1}, \ldots, Y_{m})
\end{equation}

#### Computing p-values

* Give exercise, compute p-values for Wilcoxon test where
we have two populations. both are Normally distributed
with mean zero but different variances.


## One Sample Tests

### The Sign Test

* Suppose we have observations $W_{1}, \ldots, W_{n}$ which arise from the following model
\begin{equation}
W_{i} = \theta + \varepsilon_{i}, \nonumber 
\end{equation}
where $\varepsilon_{i}$ are iid random variables each with distribution function $F$
that is assumed to have a median of zero.



### The Signed-Rank Wilcoxon Test

## Comparisons with Parametric Tests

## Thinking about Rank statistics more generally




