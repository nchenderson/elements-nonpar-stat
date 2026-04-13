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

