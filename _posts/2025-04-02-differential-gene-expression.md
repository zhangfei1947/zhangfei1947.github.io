---
title: Differential Gene Expression
author: me
date: 2025-04-02 09:59:00 -0500
categories: [Statistics, Algorithm]
tags: [statistics, algorithm, RNA-Seq]
math: true
render_with_liquid: false
---

## Binomial Distribution and Poission Distribution
1. Bernoulli trail success probability $p$
1. $n$ times independent trails
1. obtaining $k$ successes

$$
P(X = K) = \binom{n}{k}p^k(1-p)^{n-k},\ k = 0,1,2,...,n
$$

Mean: $\mu = np$

$$
E[X] = E[\displaystyle\sum_{i=1}^{n} X_i ] = \displaystyle\sum_{i=1}^{n} E[X_i] = np
$$

Variance: $Var(X) = np(1-p)$

\begin{equation}
\begin{split}
	Var(X) & = E[X^2] - (E[X])^2 \\
	& = E[X(X-1)] + E[X] - (E[X])^2 \\
	& = \displaystyle\sum_{k=0}^{n} k(k-1)\binom{n}{k}p^kq^{(n-k)} + np - (np)^2 \\
	& = k(k-1)p^2\displaystyle\sum_{k=0}^{n} \binom{n}{k}p^{k-2}q^{(n-k)} + np - (np)^2 \\
\end{split}
\end{equation}

## Negative Binomial Distribution
1. Bernoulli trail success probability $p$
1. $X$ times failures before the $r$-th success

$$
P(X=x) = \binom{k+r-1}{r-1}p^r(1-p)^k
$$

Mean: $\mu = \frac{r(1-p)}{p}$

Variance: $Var(X) = \frac{r(1-p)}{p^2}$



### Two common parameterization

## Overdispersion of gene expression

## Parameterization of NB model

$$
Var(y) = \mu + \phi*\mu^2
$$

## Estimation of Dispersion $\phi$

$$
Y_{gi} \sim NB(\mu_{gi}, \phi_g)
$$

## Differential Test
### DESeq2
**Wald Test**

**Likelihood Ration Test, LRT**

### edgeR
**Quasi-Likelihood F Test, QLF**

## Multiple Tests Correction
**Benjamini-Hochberg correction**

