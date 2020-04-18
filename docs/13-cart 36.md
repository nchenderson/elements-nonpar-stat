# (PART) Nonparametric Regression: Part II {-} 

# Decision Trees and CART {#decision-tree}

## Introduction

* Let's think about the regressogram estimate again. 

* The regressogram estimate of the regression function is a piecewise constant function that is constant within each of $p$ "bins"
\begin{equation}
\hat{m}_{h_{n}}^{R}(x) = \frac{1}{a}\sum_{i=1}^{n} Y_{i} I(x_{i} \in B_{k}), \qquad \textrm{ if } x \in B_{k} \nonumber
\end{equation}

* Figure \@ref(fig:cart-motivate) shows an example of a regressogram estimate with 4 bins. 

<div class="figure">
<img src="13-cart_files/figure-html/cart-motivate-1.png" alt="Regressogram estimate with 4 bins." width="672" />
<p class="caption">(\#fig:cart-motivate)Regressogram estimate with 4 bins.</p>
</div>

* For the data shown in \@ref(fig:cart-motivate), 


