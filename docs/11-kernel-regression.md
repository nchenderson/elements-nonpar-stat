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

---

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

---

* The estimate $\hat{m}_{h_{n}}^{R}(x)$ of the regression function is called the **regressogram**.

* The intuition for this estimate is: if $x \in B_{k}$,
then taking an average of the reponses for $x_{i}$ in a small bin containing $x$ 
should give us a reasonable approximation for the expectation of $Y_{i}$ given that $x_{i} = x$.

* Another way of looking at the regressogram is to note that for $x \in B_{k}$
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

---

<div class="figure">
<img src="11-kernel-regression_files/figure-html/unnamed-chunk-1-1.png" alt="Framingham Data. Regressogram estimate for a regression model with diastolic blood pressure as the response and age as the covariate. Ages from 31-71 were separated into bins of width 5 years." width="672" />
<p class="caption">(\#fig:unnamed-chunk-1)Framingham Data. Regressogram estimate for a regression model with diastolic blood pressure as the response and age as the covariate. Ages from 31-71 were separated into bins of width 5 years.</p>
</div>

---

* **Exercise 11.1** Let $\hat{\mathbf{Y}} = (\hat{Y}_{1}, \ldots, \hat{Y}_{n})$ denote
the vector of "fitted values" from a regressogram estimate that has $D_{n}$ bins. 
That is, $\hat{Y}_{i}$ is defined as
\begin{equation}
\hat{Y}_{i} = \hat{m}_{h_{n}}^{R}(x_{i}) \nonumber
\end{equation}
If $\mathbf{Y} = (Y_{1}, \ldots, Y_{n})$, show that you can express $\hat{\mathbf{Y}}$ as 
\begin{equation}
\hat{\mathbf{Y}} = \mathbf{A}\mathbf{Y}, \nonumber
\end{equation}
for an appropriately chosen $n \times n$ matrix $\mathbf{A}$.
What is the value of $\textrm{tr}(\mathbf{A})$?

---


### The Local Average Estimator

* The regressogram can be thought of as a regression analogue of the histogram.

* The local average estimator can be thought of as a regression analogue of the
"box-type" density estimator that we described in Chapter 8.


## Additional Reading

* Additional reading which covers the material discussed in this chapter includes:
    + Chapter 4 from @hardle2012
    + Chapter 5 from @wasserman2006
