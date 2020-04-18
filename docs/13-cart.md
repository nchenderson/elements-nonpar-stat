# (PART) Nonparametric Regression: Part II {-} 

# Decision Trees and CART {#decision-tree}

## Introduction

* Let's think about the regressogram estimate again. 

* The regressogram estimate of the regression function is a piecewise constant function that is constant within each of $p$ "bins"
\begin{equation}
\hat{m}_{h_{n}}^{R}(x) = \frac{1}{a}\sum_{i=1}^{n} Y_{i} I(x_{i} \in B_{k}), \qquad \textrm{ if } x \in B_{k} \nonumber
\end{equation}

* Figure \@ref(fig:cart-motivate) shows an example of a regressogram estimate with 3 bins. 

<div class="figure">
<img src="13-cart_files/figure-html/cart-motivate-1.png" alt="Regressogram estimate with the 3 bins: [0,1/3), [1/3, 2/3), [2/3, 1)." width="672" />
<p class="caption">(\#fig:cart-motivate)Regressogram estimate with the 3 bins: [0,1/3), [1/3, 2/3), [2/3, 1).</p>
</div>

---

* Suppose we were forced to combined two adjacent bins to estimate the regression function. For the data shown in \@ref(fig:cart-motivate), 
which two bins should we combine if we were forced to do so? The only options here are to combine bins 1 and 2 or to combine bins 2 and 3.

* I would say we should combine bins 1 and 2. Look at Figures \@ref(fig:cart-motivate2) and \@ref(fig:cart-motivate3) for a comparison of these two choices.

* The responses $Y_{i}$ change much more over the range of the third bin than they do over the first and second bins. Hence, an 
intercept model for the first two bins is not all that bad. 

* In contrast, an intercept model for the last two bins is a terrible model. 

---

* This intuition for why the choice of bins $[0,2/3), [2/3, 1)$ is better than the choice of bins $[0,1/3), [1/3, 1)$ can 
be formalized by considering the within-bin variance of $Y_{i}$.

* If we have $p$ bins, the within-bin variance of $Y_{i}$ would be
\begin{equation}
WBVar = \sum_{k=1}^{p} \frac{1}{n_{k}-1}\sum_{i=1} (Y_{i} - \bar{Y}_{k})^{2}I(x_{i} \in B_{k}) \nonumber
\end{equation}


<div class="figure">
<img src="13-cart_files/figure-html/cart-motivate2-1.png" alt="Regressogram estimate with the 2 bins: [0,1/3), [1/3, 1)." width="672" />
<p class="caption">(\#fig:cart-motivate2)Regressogram estimate with the 2 bins: [0,1/3), [1/3, 1).</p>
</div>

<div class="figure">
<img src="13-cart_files/figure-html/cart-motivate3-1.png" alt="Regressogram estimate with the 2 bins: [0,2/3), [2/3, 1)." width="672" />
<p class="caption">(\#fig:cart-motivate3)Regressogram estimate with the 2 bins: [0,2/3), [2/3, 1).</p>
</div>



```r
var(yy[ind1 | ind2])
```

```
## [1] 0.07720792
```

```r
var(yy[ind3])
```

```
## [1] 0.4485642
```

```r
var(yy[ind2 | ind3])
```

```
## [1] 0.9214935
```

```r
var(yy[ind1])
```

```
## [1] 0.05972719
```










