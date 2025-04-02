---
title: Differential gene expression
author: me
date: 2025-04-02 09:59:00 -0500
categories: [Statistics, Algorithm]
tags: [statistics, algorithm, RNA-Seq]
math: true
render_with_liquid: false
---

## Binomial Distribution
1. Bernoulli trail success probability: $p$
1. $n$ times independent trails
1. obtaining $k$ successes
$$
P(X = K) = \binom{n}{k}p^k(1-p)^{n-k}
$$

Mean: $\mu = np$

Variance: $\sigma^2 = npq$

## Negative Binomial Distribution
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

