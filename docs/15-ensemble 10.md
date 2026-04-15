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
## [1]	train-rmse:130.780291±0.582002	test-rmse:130.788142±2.388444 
## [51]	train-rmse:84.854008±0.403340	test-rmse:85.275384±1.595941 
## [101]	train-rmse:74.971515±0.277799	test-rmse:75.524741±1.094459 
## [151]	train-rmse:71.671575±0.485694	test-rmse:72.507682±0.921721 
## [201]	train-rmse:68.839892±0.446413	test-rmse:69.800556±1.065498 
## [251]	train-rmse:66.963632±0.371290	test-rmse:68.039514±1.155026 
## [301]	train-rmse:65.578527±0.434747	test-rmse:66.764971±1.202784 
## [351]	train-rmse:64.327758±0.528381	test-rmse:65.586123±1.146745 
## [401]	train-rmse:63.019676±0.460593	test-rmse:64.383894±1.369827 
## [451]	train-rmse:61.690550±0.496092	test-rmse:63.140779±1.493872 
## [501]	train-rmse:60.506610±0.448749	test-rmse:62.030273±1.657523 
## [551]	train-rmse:59.407744±0.602749	test-rmse:60.992809±1.532014 
## [601]	train-rmse:58.490998±0.499451	test-rmse:60.142795±1.543443 
## [651]	train-rmse:57.726195±0.513933	test-rmse:59.473072±1.559282 
## [701]	train-rmse:57.120484±0.504649	test-rmse:58.945631±1.566263 
## [751]	train-rmse:56.604565±0.490694	test-rmse:58.476718±1.499277 
## [801]	train-rmse:56.068782±0.530417	test-rmse:57.968260±1.441840 
## [851]	train-rmse:55.615336±0.508406	test-rmse:57.563509±1.436265 
## [901]	train-rmse:55.200157±0.450105	test-rmse:57.214875±1.439893 
## [951]	train-rmse:54.814117±0.444433	test-rmse:56.864085±1.411996 
## [1000]	train-rmse:54.468861±0.440165	test-rmse:56.586989±1.401586
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
## 1:     1        130.7803      0.5820018       130.7881      2.388444
## 2:     2        128.0044      0.5740515       128.0144      2.412024
## 3:     3        125.4426      0.5673276       125.4800      2.427631
## 4:     4        123.0701      0.5636404       123.1289      2.447933
## 5:     5        120.8712      0.5617252       120.9460      2.448483
## 6:     6        118.8365      0.5567243       118.9606      2.421264
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
