---
title: Algorithm of JTK_CYCLE
author: me
date: 2025-04-01 09:59:00 -0500
categories: [Statistics, Algorithm]
tags: [statistics, algorithm, circadian]
math: true
render_with_liquid: false
---

Hughes, Michael E., John B. Hogenesch, and Karl Kornacker. "**JTK_CYCLE: an efficient nonparametric algorithm for detecting rhythmic components in genome-scale data sets.**" _Journal of biological rhythms_ 25.5 (2010): 372-380.

doi: [10.1177/0748730410379711](https://doi.org/10.1177/0748730410379711)

JTK_CYCLE is a classic and widely used software to identify rhythmic expressed genes through high throughput technologies.

## 1. Jonckheere-Terpstra (JT) Test


$$
U = \sum_{i=1}^{k-1} \sum_{j=i+1}^k U_{ij}
$$

For large sample, under $H_0$

$$
Z = \frac{U - E(U)}{\sqrt{\text{Var}(U)}}
$$

## 2. Kendall's $\tau$


## 3. JTK


## 4. Harding Algorithm

