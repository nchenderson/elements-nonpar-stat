# Tree Ensemble Methods for Prediction {#ensemble}
 

## Motivation

* A major issue with single-decision tree models is that they can have **high variance**.
    + At least among decision trees that are allowed to be flexible enough to have modest bias.
    
* Small perturbations in the data can often have substantial impacts on the final estimated tree.

* This can occur, for example, because a small perturbation can impact the first variable chosen to split on. 
    + This early split further influences splitting variables and splitting points chosen further down the tree.
    
* For single-decision trees, you can have substantial changes in the fitted values when moving from
one region to another.

* You can often have small number of observations in some terminal nodes.

---

* Using **ensembles of trees** has been found 

* Two well-known tree ensemble methods:
    + **Random Forests**
    + **Gradient Boosting** with Trees

## Gradient Boosting with Regression Trees

* Boosting is an ensemble method that takes a **sequential appraoch** to building an ensemble of trees.

* Uses a somewhat similar idea to **forward selection** in linear regression.

* Each tree is typically a "shallow tree" (weak learner).
    + Each tree, by itself, often **does not** have strong performance.
    + Each tree: high bias, but low variance.
    + Aggregating all trees leads to stronger performance.

- Idea: Start with a **single tree**
    -   Add more trees one at a time.
    -   Each new tree tries to improve on the previous ones.

---

* **Goal:** We want to estimate a regression function $f_{M}(\mathbf{x})$ that can be expressed as the sum of M trees: 
\begin{equation}
    f_{M}(\mathbf{x}) = \sum_{m=1}^{M} T(\mathbf{x}; \Theta_{m}) 
    = \sum_{m=1}^{M} \sum_{h=1}^{H} \mu_{hm} I( \mathbf{x} \in A_{hm} ). \nonumber 
\end{equation}

* To describe this procedure, let $T(\mathbf{x}; \Theta_{m})$ denote the "fitted value" for the input covariate vector $\mathbf{x}$ when using a tree with parameters $\Theta_{m}$.

* Parameters of the $m^{th}$ tree: $\Theta_{m} = \{ A_{hm}, \mu_{hm} \}_{h=1}^{H}$, where
    + $H$ disjoint \`\`terminal regions'' $A_{1m}, \ldots, A_{Hm}$
    + $\mu_{hm}$ output assigned to region $A_{hm}$, i.e., $T(\mathbf{x}; \Theta_{m}) = \mu_{hm}$ whenever $\mathbf{x} \in A_{hm}$.

## Estimating the Tree in each Iteration

* Strategy: estimate parameters $\Theta_{1}, \ldots, \Theta_{M}$ in a **forward, stagewise manner**.

* First, estimate $\widehat{\Theta}_{1}$ by targeting the **loss function** $L(\mathbf{y}, T(\mathbf{X}; \Theta_{1}))$

* For example, the squared-error loss would be 
\begin{equation}
L(\mathbf{y}, T(\mathbf{X}; \Theta_{1})) = \sum_{i=1}^{n} \Big( y_{i} - T(\mathbf{x}_{i}; \Theta_{1}))
\end{equation}

* Remaining tree parameters are estimated by looking at 
\begin{equation}
    L_{I}\Big( \mathbf{Y}, \sum_{k=1}^{m-1} T(\mathbf{X}; \widehat{\Theta}_{k})  + T(\mathbf{X}; \Theta_{m}) \Big).
    \label{eq:direct_boost_loss}
\end{equation}

* $m^{th}$ iteration of boosting: estimate $\Theta_{m}$

---

* Estimation of $\Theta_{m}$ is performed in **two stages**

* **Stage 1**: Calculate **terminal regions** $A_{hm}$ by fitting the (generalized) residuals $r_{i,m-1}$ \begin{equation}
    r_{i,m-1} = -\Big[ \frac{\partial L(y_{i}, f(\mathbf{x}_{i}) }{ \partial f(\mathbf{x}_{i})} \Big]_{f = f_{m-1}}
    \end{equation}

* Fit a regression tree using the residuals $r_{1,m-1}, \ldots, r_{n,m-1}$ as outcomes.

* The number of terminal regions is often fixed to a small number or chosen to have a maximum number of terminal regions (such as 6).

---

* **Stage 2**: Estimate terminal values $\mu_{hm}$ within each region $A_{hm}$, by minimizing:
\begin{equation}
    \mu_{hm} = \arg\min_{\alpha} \sum_{x_{i} \in A_{hm}} L\Big(y_{i}, f_{m-1}(\mathbf{x_{i}}) + \alpha \Big)
\end{equation}

### Adding Trees and Learning Rate

* Updating the boosting model by using **shrinkage** is common in most implementations of boosting.

* Obtain $f_{m}(\mathbf{x})$ from $f_{m-1}(\mathbf{x})$ with \begin{equation}
    f_{m}(\mathbf{x}) = f_{m-1}(\mathbf{x}) + \eta T(\mathbf{x}, \widehat{\Theta}_{m}),
    \label{eq:shrink_boost_update}
    \end{equation}

* The shrinkage term $\eta \in (0, 1)$ is often referred to as the **\`\`learning rate''**

-   Choosing a relatively small value of $\gamma$ typically leads to improved out-of-sample predictive performance.
    + Usually, $\eta \leq 0.3$.

---

* Boosting often requires tuning of the hyperparameters.

* Main hyperparameters that require tuning:
    + Learning rate $\eta$.
    + Number of trees $M$.

### Choice of Loss Function

* An advantage of gradient boosting is that it can work for a general loss function (as long as the loss is differentiable with respect to $f$).

* Because of this, gradient boosting can be applied to many types of data.

* **Binary outcomes**: 
\begin{equation}
    L(\mathbf{y}, f) = \sum_{i=1}^{n}L(y_{i}, f(\mathbf{x}_{i})) = -\sum_{i=1}^{n}y_{i}\log( f(\mathbf{x}_{i})) - \sum_{i=1}^{n}(1 - y_{i})\log(1 - f(\mathbf{x}_{i}))
\end{equation}

## XGBoost

* `xgboost` is a particular version of gradient boosting that is one of the more popular boosting procedures.


``` r
library(xgboost)
```

### XGBoost Example:

* To try `xgboost`, let's use the `Bikeshare` data again




* Use the function `xgb.DMatrix` to store in the standard `xgboost` format


``` r
X <- model.matrix(bikers ~ . - casual - registered, data=Bikeshare)
dtrain <- xgb.DMatrix(data = X, label = Bikeshare$bikers)
```


* Setup the **loss function** and tuning parameters


``` r
params <- list(
  objective = "reg:squarederror", 
  eta = 0.05,                   
  max_depth = 2                
)
```

* Use these tuning parameters to evaluate cross-validation loss for up to 1000 trees


``` r
cv_results <- xgb.cv(params = params,
  data = dtrain,
  nrounds = 1000,               # Set an artificially high ceiling
  nfold = 5,
  early_stopping_rounds = 20,   # Stop if validation loss doesn't improve for 20 trees
  print_every_n = 50,           # Print progress every 50 trees
  verbose = TRUE)
```

```
## Multiple eval metrics are present. Will use test_rmse for early stopping.
## Will train until test_rmse hasn't improved in 20 rounds.
## 
## [1]	train-rmse:130.778523±0.719314	test-rmse:130.796331±2.954253 
## [51]	train-rmse:84.659287±0.474192	test-rmse:84.987781±2.136134 
## [101]	train-rmse:74.981565±0.453892	test-rmse:75.427287±1.822934 
## [151]	train-rmse:71.548004±0.286702	test-rmse:72.224634±1.773907 
## [201]	train-rmse:68.743161±0.345865	test-rmse:69.621760±1.331848 
## [251]	train-rmse:66.913853±0.459587	test-rmse:67.924167±1.094197 
## [301]	train-rmse:65.554391±0.491843	test-rmse:66.689281±1.116398 
## [351]	train-rmse:64.262514±0.498881	test-rmse:65.495011±1.268142 
## [401]	train-rmse:63.069767±0.332063	test-rmse:64.418160±1.322234 
## [451]	train-rmse:61.627847±0.390297	test-rmse:63.064030±1.438689 
## [501]	train-rmse:60.362758±0.458443	test-rmse:61.895606±1.381340 
## [551]	train-rmse:59.397808±0.532251	test-rmse:61.005388±1.252942 
## [601]	train-rmse:58.504221±0.624739	test-rmse:60.185085±1.303737 
## [651]	train-rmse:57.734228±0.564315	test-rmse:59.499109±1.251125 
## [701]	train-rmse:57.068382±0.559492	test-rmse:58.891514±1.169635 
## [751]	train-rmse:56.550684±0.564915	test-rmse:58.436166±1.135934 
## [801]	train-rmse:56.086250±0.564035	test-rmse:58.035134±1.114464 
## [851]	train-rmse:55.700891±0.573919	test-rmse:57.722340±1.102247 
## [901]	train-rmse:55.305584±0.572521	test-rmse:57.397843±1.104590 
## [951]	train-rmse:54.966606±0.537094	test-rmse:57.139070±1.113257 
## [1000]	train-rmse:54.621998±0.531869	test-rmse:56.848063±1.143896
```

``` r
## Get best number of trees
best_ntrees <- cv_results$niter
print(best_ntrees)
```

```
## [1] 1000
```

* The `evaluation_log`, gives you the cross-validation estimate of test error.


``` r
cv_est <- cv_results$evaluation_log
print(head(cv_est))
```

```
##     iter train_rmse_mean train_rmse_std test_rmse_mean test_rmse_std
##    <int>           <num>          <num>          <num>         <num>
## 1:     1        130.7785      0.7193142       130.7963      2.954253
## 2:     2        128.0033      0.7059292       128.0246      2.950007
## 3:     3        125.4428      0.6942594       125.4741      2.956300
## 4:     4        123.0729      0.6819842       123.1256      2.953471
## 5:     5        120.8760      0.6672395       120.9262      2.926570
## 6:     6        118.8431      0.6536342       118.9114      2.940425
```

* Let's plot the estimate of test error versus the number of trees


``` r
plot(cv_est$iter, cv_est$test_cox_nloglik_mean, xlab="Number of Trees",
     ylab="Estimate of Test Error", lwd=2, cex.lab=1.8, cex.axis=1.8)
lines(cv_est$iter, cv_est$test_cox_nloglik_mean, lwd=2)
```

* `xgboost` **final model:**

``` r
fmodel <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = best_ntrees, ## Use best number of trees found.
  verbose = FALSE 
)
```

* For variable importance, you can use the `xgb.importance` function:

* Focus first on the `Gain` measure


``` r
importance_matrix <- xgb.importance(model = fmodel)
print(importance_matrix)
```
