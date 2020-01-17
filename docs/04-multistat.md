# Rank Tests for Multiple Groups {#krusk-wallis}

* We can roughly think of the tests discussed in Chapter 3
as being related to the parametric tests shown in the table below.


* The **Kruskal-Wallis** test can be though of as the
nonparametric analogue of one-way analysis of variance (ANOVA).

* For $K$ groups, one-way ANOVA considers the analysis of data
arising from the following model
\begin{equation}
Y_{kj} = \mu_{k} + \varepsilon_{kj}, \qquad j=1,\ldots, n_{k}; k=1,\ldots,K  \nonumber
\end{equation}
where it is often assumed that $\varepsilon_{kj} \sim \textrm{Normal}(0, \sigma^{2})$.


* Usually, the one-way ANOVA hypothesis of interest is something like
\begin{equation}
H_{0}: \mu_{1} = \mu_{2} = \ldots = \mu_{K}
\end{equation}

* A test of the hypothesis () is based on decomposing the observed variation in
the responses
\begin{equation}
\sum_{k}\sum_{j} (Y_{kj} - \bar{Y})^{2} = \sum_{k} \sum_{j} (\bar{Y}_{k.} - \bar{Y}_{..})^{2}
+ \sum_{k}\sum_{j}(\bar{Y}_{kj} - \bar{Y}_{k.})^{2}
\end{equation}

* Large values 

---

* Instead of assuming () for the data $Y_{kj}$, nonparametric way of thinking
about this problem is to instead only assume that
\begin{equation}
Y_{kj} \sim F_{k}
\end{equation}
That is, $Y_{k1}, Y_{k2}, \ldots, Y_{kn_{k}}$ is an i.i.d. sample from $F_{k}$. 


* A nonparametric version of the one-way ANOVA hypothesis is that
\begin{equation}
F_{1} = F_{2} = \ldots = F_{K}
\end{equation}

---

* The Kruskall-Wallis test statistic is defined as
\begin{equation}
KW_{n} = \frac{12}{N(N-1)}\sum_{k=1}^{K} n_{k}\Big( \bar{R}_{k.} - \frac{N + 1}{2} \Big)
\end{equation}
