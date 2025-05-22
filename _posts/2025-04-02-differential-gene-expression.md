---
title: Differential Gene Expression
author: me
date: 2025-04-02 09:59:00 -0500
categories: [Statistics, Algorithm]
tags: [statistics, algorithm, RNA-Seq]
math: true
render_with_liquid: false
---
## Step 1. Normalization of Library Size

![Sequencing depth](/assets/img/DGE/SeqDepthVaries.png){: w="350" }

### 1. Median of Ratios (DESeq2)

1. Calculate geometric mean of all genes across all samples as a "pseudo-reference"
1. Calculate ratios of gene readcount to it's geometric mean for each sample
1. Let the median of ratios as size factor for each sample

$$ s_i = \text{median}_g \left( \frac{K_{gi}}{\left( \prod_{j} K_{gj} \right)^{1/m}} \right) $$

![Norm seq depth](/assets/img/DGE/MoR_Norm.png){: w="350" }

### 2. Trimmed Mean of M-values (edgeR)




## Step 2. Modeling Feature's Read Count

### 1. Binomial Distribution
1. Bernoulli trail success probability $p$
1. $n$ times independent trails
1. Probability of obtaining $X$ successes

$$
P(X = x) = \binom{n}{k}p^x(1-p)^{n-x}\ \ x = 0,1,2,...,n
$$

**Mean:** $\mu = np$

$$
E[X] = E[\displaystyle\sum_{i=1}^{n} X_i ] = \displaystyle\sum_{i=1}^{n} E[X_i] = np
$$

**Variance:** $Var(X) = np(1-p)$

$$
\begin{align*}
	Var(X) & = E[X^2] - (E[X])^2 \\\\
	& = E[X(X-1)] + E[X] - (E[X])^2 \\\\
	& = \displaystyle\sum_{x=0}^{n} x(x-1)\binom{n}{x}p^xq^{n-x} + np - (np)^2 \\\\
	& = \displaystyle\sum_{x=2}^{n} x(x-1) \frac{n!}{x!(n-x)!} p^xq^{n-x} + np - (np)^2 \\\\
	& = n(n-1)p^2 \displaystyle\sum_{x=2}^{n} \frac{(n-2)!}{(x-2)!(n-x)!} p^{x-2}q^{n-x} + np - (np)^2 \\\\
	& = (np)^2 - np^2 + np - (np)^2 \\\\
	& = np(1-p)
\end{align*}
$$

### 2. Poisson Distribution
From Binomial Distribution:

1. Let $n \rightarrow \infty$
1. Let $p \rightarrow 0$
1. Meanwhile, let the expectation $\mu = np$ be a constant $\lambda$

$$
P(X=x) = \frac{e^{-\lambda}\lambda^x}{x!}
$$

**Mean:** $\mu = \lambda$

**Variance:** $Var(X) = \lambda$


### 3. Negative Binomial Distribution
From Binomial Distribution:

1. Probability of get $X$ times failures before the $r$-th success

$$
P(X=x) = \binom{x+r-1}{r-1}p^r(1-p)^x\ \ x=0,1,2... \\\\
\binom{x+r-1}{r-1} = (-1)^x \binom{-r}{x}
$$

**Mean:** $\mu = \frac{r(1-p)}{p}$

**Variance:** $Var(X) = \frac{r(1-p)}{p^2}$


#### 3.1. NB is a mixture model of Poisson and Gamma

$$
X|\lambda \sim Poisson(\lambda): \ \ P(X=x) = \frac{e^{-\lambda}\lambda^x}{x!}
$$

$$
\lambda \sim Gamma(\alpha, \beta): \ \ f(\lambda) = \frac{\beta^\alpha}{\Gamma(\alpha)}\lambda^{\alpha-1}e^{-\beta\lambda}
$$

$$
P(X=x) = \int_0^\infty P(X=x|\lambda) f(\lambda;\alpha,\beta) \, d\lambda = ... = \binom{x+\alpha-1}{x}p^x(1-p)x \ \ \ \ (p=\frac{\beta}{1+\beta})
$$

#### 3.2. Other parameterizations of NB

1. Let $X$ be total tail count when get $r$-th success

$$
P(X=x|r,p) = \binom{x-1}{r-1}p^r(1-p)^{x-r}\ \ \ x=r, r+1, ...
$$

2. Notice that $Var(X)=E(X)/p$, it can be used to modeling overdispersion data. Let expectation $E(X) = \mu$

$$
\mu = \frac{r(1-p)}{p} \\\\
p = \frac{r}{\mu+r} \\\\
1-p = \frac{\mu}{\mu+r} \\\\
\begin{align*}
Var(X) &= \frac{r(1-p)}{p^2} \\\\
&= \mu + \frac{\mu^2}{\alpha} \ \ \ (let\ \alpha=r) \\\\
&= \mu + \phi*\mu^2
\end{align*}
$$

$\alpha$: dispersion coefficient, when $\alpha \downarrow$, dispersion level $\uparrow$


### Overdispersion of RNA-Seq gene readcount

![Mean-Variance plot](/assets/img/DGE/Overdispersion.png){: w="400" }

### Negative binomial model applied

$$
K_{gi} \sim NB(\mu_{gi}, \alpha_g)
$$

$K_{gi}$: read count of gene $g$ in sample $i$

$\mu_{gi}$: expected read count of gene $g$ in sample $i$

$\alpha_g$: dispersion parameter, indicate variation level of gene $g$

Before performing statistical tests between $\mu_{gi}$, obtaining a better estimation of $\alpha_g$ ensures more accurate NB modeling.

### Estimation of Dispersion$\phi$ (DESeq2)

1. Initial (Per-Gene) Estimation: calculates a Maximum Likelihood Estimate (MLE) of the dispersion for each individual gene

![Init MLE](/assets/img/DGE/InitialMLE.png){: w="400" }

$$
K_{gi} \sim NB(\mu_{gi}, \alpha_g) \\\\
P(X=x) = \binom{x+r-1}{r-1}p^r(1-p)^x \\\\
x = k \\\\
r = \frac{1}{\alpha} \\\\
p = \frac{1/\alpha}{\mu_{gi}+1/\alpha} \\\\
$$

$$
\begin{align*}
P(K_{gi} = k) &= \binom{k + \frac{1}{\alpha_g} - 1}{k} \left( \frac{1}{1 + \alpha_g \mu_{gi}} \right)^{\frac{1}{\alpha_g}} \left( \frac{\alpha_g \mu_{gi}}{1 + \alpha_g \mu_{gi}} \right)^k \\\\
&= \frac{\Gamma(k + \frac{1}{\alpha_g})}{\Gamma(k+1)\Gamma(\frac{1}{\alpha_g})} \left( \frac{1}{1 + \alpha_g \mu_{gi}} \right)^{\frac{1}{\alpha_g}} \left( \frac{\alpha_g \mu_{gi}}{1 + \alpha_g \mu_{gi}} \right)^k
\end{align*}
$$

Likelihood function:

$$
L(\alpha_g) = \prod_{i=1}^{n} P(K_{gi}) \\\\
\ell(\alpha_g) = \sum_{i=1}^{n} log P(K_{gi}) \\\\
\ell(\alpha_g) = \sum_{i=1}^{n} \left[ \log \Gamma\left(K_{gi} + \frac{1}{\alpha_g} \right) - \log \Gamma\left( \frac{1}{\alpha_g} \right) - K_{gi} \log(1 + \alpha_g \mu_{gi}) + K_{gi} \log \alpha_g \mu_{gi} \right]
$$


1. Fitting a Global Trend Curve: fits a curve between the mean expression level of a gene and its dispersion

![Trend Curve](/assets/img/DGE/TrendCurve.png){: w="400" }
LOESS（Locally Weighted Scatterplot Smoothing)

$$
\log_2(\alpha_g) = \text{LOESS}(\log_2(\mu_g))
$$

1. Shrinkage using Empirical Bayes: initial per-gene dispersion estimates are "shrunk" towards the fitted global trend curve

![Shrink](/assets/img/DGE/Shrink2Curve.png){: w="400" }

$$
\hat{\alpha}_g^{\text{(shrunk)}} = w_g \cdot \hat{\alpha}_g^{(MLE)} + (1 - w_g) \cdot \alpha_g^{\text{(trend)}}
$$


## Step 3. Statistical Tests
### DESeq2
**Wald Test**

**Likelihood Ration Test, LRT**

### edgeR
**Quasi-Likelihood F Test, QLF**

### Multiple Tests Correction
**Benjamini-Hochberg correction**


