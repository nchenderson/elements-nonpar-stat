# Permutation Tests {#permutation}

## Notation

* A **permutation** $\pi$ of a set $S$ is a function $\pi: S \longrightarrow S$ is a
function that is both one-to-one and onto.

* We will usually think of $S$ as the set of observation indices in which case
$S = \{1, \ldots, N\}$ for sample size $N$. 

* Each permutation $\pi$ of $S = \{1, \ldots, N\}$ defines a particular ordering of the elements of $S$.
For this reason, a permutation is often expressed as the following ordered list
\begin{equation}
\pi = \big( \pi(1), \pi(2), \ldots, \pi(N)  \big)
\end{equation}

* In other words, we can think of a permutation of $S$
as a particular ordering of the elements of $S$.


* For example, if $S = \{1,2,3\}$, and $\pi_{1}$ is a permuation of $S$
defined as $\pi_{1}(1) = 3$, $\pi_{1}(2) = 1$, $\pi_{1}(3) = (2)$, then
this permutation expressed as an ordered list would be
\begin{equation}
\pi_{1} = (3, 1, 2)
\end{equation}

* There are $5$ other possible permutations of $S$: 
\begin{eqnarray}
\pi_{2} &=& (1,2,3) \nonumber \\
\pi_{3} &=& (2,1,3) \nonumber \\
\pi_{4} &=& (1,3,2) \nonumber \\
\pi_{5} &=& (3,2,1) \nonumber \\
\pi_{6} &=& (2, 3, 1) \nonumber  
\end{eqnarray}

* If $S$ has $N$ distinct elements, there are $N!$ possible permutations of $S$. 

## Permutation Tests for the Two-Sample Problem

* A permutation test is motivated by the following reasoning.
    + If there is no real difference between the two groups, 
      there is noting "special" about the difference in means
      between the two groups.
    + The difference 
      
      

