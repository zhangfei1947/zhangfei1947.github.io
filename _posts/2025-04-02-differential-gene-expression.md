---
title: Differential Gene Expression
author: me
date: 2025-04-02 09:59:00 -0500
categories: [Statistics, Algorithm]
tags: [statistics, algorithm, RNA-Seq]
math: true
render_with_liquid: false
---

## Binomial Distribution
1. Bernoulli trail success probability $p$
1. $n$ times independent trails
1. Probability of obtaining $X$ successes

$$
P(X = x) = \binom{n}{k}p^x(1-p)^{n-x},\ \ x = 0,1,2,...,n
$$

**Mean:** $\mu = np$

$$
E[X] = E[\displaystyle\sum_{i=1}^{n} X_i ] = \displaystyle\sum_{i=1}^{n} E[X_i] = np
$$

**Variance:** $Var(X) = np(1-p)$

\begin{equation}
\begin{split}
	Var(X) & = E[X^2] - (E[X])^2 \\
	& = E[X(X-1)] + E[X] - (E[X])^2 \\
	& = \displaystyle\sum_{x=0}^{n} x(x-1)\binom{n}{x}p^xq^{n-x} + np - (np)^2 \\
	& = \displaystyle\sum_{x=2}^{n} x(x-1) \frac{n!}{x!(n-x)!} p^xq^{n-x} + np - (np)^2 \\
	& = n(n-1)p^2 \displaystyle\sum_{x=2}^{n} \frac{(n-2)!}{(x-2)!(n-x)!} p^{x-2}q^{n-x} + np - (np)^2 \\
	& = (np)^2 - np^2 + np - (np)^2 \\
	& = np(1-p)
\end{split}
\end{equation}

## Poisson Distribution
From Binomial Distribution:

1. Let n &rarr; &infin;
2. Let p &rarr; 0
3. Meanwhile, let the expectation $\mu = np$ be a constant $np = \lambda$

$$
P(X=x) = \frac{e^{-\lambda}\lambda^k}{k!}
$$

**Mean:** $\mu = \lambda$

**Variance:** $Var(X) = \lambda$


## Negative Binomial Distribution
1. Bernoulli trail success probability $p$
1. Probability of get $X$ times failures before the $r$-th success

$$
P(X=x) = \binom{x+r-1}{r-1}p^r(1-p)^x,\ \ x=0,1,2... \\
\binom{x+r-1}{r-1} = (-1)^x \binom{-r}{x}
$$

**Mean:** $\mu = \frac{r(1-p)}{p}$

**Variance:** $Var(X) = \frac{r(1-p)}{p^2}$



### Two common parameterization
Let $X$ be total tail count when get $r$-th success



## Overdispersion of RNA-Seq gene readcount

![Mean-Variance plot](/assets/img/Overdispersion.png)

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


$$s_i = \text{median}_g \left( \frac{K_{gi}}{\left( \prod_{j} K_{gj} \right)^{1/m}} } \right)$$


