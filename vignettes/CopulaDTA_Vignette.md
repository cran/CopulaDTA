---
title: "CopulaDTA Package Vignette"
author: "Victoria N Nyaga"
date: "2017-10-11"
output: rmarkdown::html_vignette
bibliography: References.bib
vignette: >
  %\VignetteIndexEntry{CopulaDTA Package Vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The current statistical procedures implemented in statistical software packages for pooling of diagnostic test accuracy data include hSROC [@hsroc] regression and the bivariate random-effects meta-analysis model (BRMA) ([@Reitsma], [@Arends], [@Chu], [@Rileya]). However, these models do not report the overall mean but rather the mean for a central study with random-effect equal to zero and have difficulties estimating the correlation between sensitivity and specificity when the number of studies in the meta-analysis is small and/or when the between-study variance is relatively large [@Rileyb].

This tutorial on advanced statistical methods for meta-analysis of diagnostic accuracy studies discusses and demonstrates Bayesian modeling using CopulaDTA [@copuladta] package in R [@R] to fit different models to obtain the meta-analytic parameter estimates. The focus is on the joint modelling of sensitivity and specificity using copula based bivariate beta distribution. Essentially, 

- we present the Bayesian approach which offers flexibility and ability to perform complex statistical modelling even with small data sets and 
- include covariate information, and 
- provide an easy to use code. 

The statistical methods are illustrated by re-analysing data of two published meta-analyses.

Modelling sensitivity and specificity using the bivariate beta distribution provides marginal as well as study-specific parameter estimates as opposed to using bivariate normal distribution (e.g., in BRMA) which only yields study-specific parameter estimates. Moreover, copula based models offer greater flexibility in modelling different correlation structures in contrast to the normal distribution which allows for only one correlation structure.

## Introduction

In a systematic review of diagnostic test accuracy, the statistical analysis section aims at estimating the average (across studies) sensitivity and specificity of a test and the variability thereof, among other measures. There tends to be a negative correlation between sensitivity and specificity, which postulates the need for correlated data models. The analysis is statistically challenging because the user 

- deals with two summary statistics, 
- has to account for correlation between sensitivity and specificity, 
- has to account for heterogeneity in sensitivity and specificity across the studies and 
- should be allowed to incorporate covariates. 

Currently, the  HSROC [@hsroc] regression or the bivariate random-effects meta-analysis model (BRMA) ([@Reitsma], [@Arends], [@Chu], [@Rileya]) are recommended for pooling of diagnostic test accuracy data. These models fit a bivariate normal distribution which allows for only one correlation structure to the logit transformed sensitivity and specificity. The resulting distribution has no closed form and therefore the mean sensitivity and specificity is only estimated after numerical integration or other approximation methods, an extra step which is rarely taken.

Within the maximum likelihood estimation methods, the BRMA and HSROC models have difficulties in convergence and estimating the correlation parameter when the number of studies in the meta-analysis are small and/or when the between-study variances are relatively large [@Takwoingi].  When the correlation is close to the boundary of its parameter space, the between study variance estimates from the BRMA are upwardly biased as they are inflated to compensate for the range restriction on the correlation parameter [@Rileyb]. This occurs because the maximum likelihood estimator truncates the between-study covariance matrix on the boundary of its parameter space, and this often occurs when the within-study variation is relatively large or the number of studies is small [@Rileya].

The BRMA and HSROC assume that the transformed data is approximately normal with constant variance, however for sensitivity and specificity and proportions in generals, the mean and variance depend on the underlying probability. Therefore, any factor affecting the probability will change the mean and the variance. This implies that the in models where the predictors affect the mean but assume a constant variance will not be adequate. 

Joint modelling of study specific sensitivity and specificity using existing or copula based bivariate beta distributions overcomes the above mentioned difficulties. Since both sensitivity and specificity take values in the interval space (0, 1), it is a more natural choice to use a beta distribution to describe their distribution across studies, without the need for any transformation. The beta distribution is conjugate to the binomial distribution and therefore it is easy to integrate out the random-effects analytically giving rise to the beta-binomial marginal distributions. Moreover no further integration is needed to obtain the meta-analytically pooled sensitivity and specificity. Previously, [@Cong] fitted separate beta-binomial models to the number of true positives and the number of false positives. While the model ignores correlation between sensitivity and specificity, [@Cong] reported that the model estimates are comparable to those from the SROC model [@Moses], the predecessor of the HSROC model. 

Ignoring the correlation would have negligible influence on the meta-analysis results when the within-study variability is large relative to the between-study variability [@Rileyc]. By utilising the correlation, we allow borrowing strength across sensitivities and specificities resulting in smaller standard errors. The use of copula based mixed models within the frequentist framework for meta-analysis of diagnostic test accuracy was recently introduced by [@Nikolo] who evaluated the joint density numerically.  

This tutorial, presents and demonstrates hierarchical mixed models for meta-analysis of diagnostic accuracy studies. In the first level of the hierarchy, given sensitivity and specificity for each study, two binomial distributions are used to describe the variation in the number of true positives and true negatives among the diseased and healthy individuals, respectively. In the second level, we  model the unobserved sensitivities and specificities using a bivariate distribution. While hierarchical models are used, the focus of meta-analysis is on the pooled average across studies and rarely on a given study estimate.

The methods are demonstrated using datasets from two previously published meta-analyses: 

- on diagnostic accuracy of telomerase in urine as a tumour marker for the diagnosis of primary bladder cancer since it is a problematic dataset that has convergence issues caused by the correlation parameter being estimated to be -1 and has no covariate [@Glas] and 
- on the comparison of the sensitivity and specificity of human papillomavirus testing (using the HC2 assay) versus repeat cytology to triage women with minor cytological cervical lesions to detect underlying cervical precancer[@Arbyn]. 
The second dataset is used to demonstrate meta-regression with one covariate which can be naturally extended to include several covariates.

## Statistical methods for meta-analysis

### Definition of copula function
A bivariate copula function describes the dependence structure between two random variables. Two random variables $X_1$ and $X_2$ are joined by a copula function C if their joint cumulative distribution function can be written as
$$F(x_1, x_2) = C(F_1 (x_1), F_2 (x_2 )), -\infty \le  x_1, x_2 \le +\infty$$

According to the theorem of [@Skylar], there exists for every bivariate (multivariate in extension) distribution a copula representation C which is unique for continuous random variables. If the joint cumulative distribution function and the two marginals are known, then the copula function can be written as
$$C(u, ~v) = F(F_1^{-1} (u), ~F_2^{-1} (v)),~ 0 \le~ u, ~v ~\le~ 1$$ 

A 2-dimensional copula is in fact simply a 2-dimensional cumulative function restricted to the unit square with standard uniform marginals. A comprehensive overview of copulas and their mathematical properties can be found in [@Nelsen]. 
To obtain the joint probability density, the joint cumulative distribution $F(x_1, x_2)$ should be differentiated to yield
$$f(x_1, ~x_2) = f_1(x_1) ~f_2(x_2 ) ~c(F_1(x_1), ~F_2 (x_2))$$
where $f_1$ and $f_2$ denote the marginal density functions and c the copula density function corresponding to the copula cumulative distribution function C. Therefore from Equation~\ref{eq:3}, a bivariate probability density can be expressed using the marginal and the copula density, given that the copula function is absolutely continuous and twice differentiable.

When the functional form of the marginal and the joint densities are known, the copula density can be derived as follows
$$c(F_1(x_1), ~F_2(x_2)) = \frac{f(x_1, ~x_2)}{f_1 (x_1 )~ f_2 (x_2 )}$$								

While our interest does not lie in finding the copula function, the equations above serve to show how one can move from the copula function to the bivariate density or vice-versa, given that the marginal densities are known. The decompositions allow for constructions of other and possible better models for the variables than would be possible if we limited ourselves to only existing standard bivariate distributions. 

We finish this section by mentioning an important implication when Sklar's theorem is extended to a meta-regression setting with covariates.  According to [@Patton], it is important that the conditioning variable remains the same for both marginal distributions and the copula, as otherwise the joint distribution might not be properly defined. This implies that covariate information should be introduced in both the marginals and the association parameters of the model.

### The hierarchical model
Since there are two sources of heterogeneity in the data, the within- and between-study variability, the parameters involved in a meta-analysis of diagnostic accuracy studies vary at two levels. For each study $i$, $i = 1, ..., n$, let $Y_{i}~=~(Y_{i1},~ Y_{i2})$  denote the true positives and true negatives, $N_{i}~=~( N_{i1},~ N_{i2})$ the diseased and healthy individuals respectively, and $\pi_{i}~ =~ (\pi_{i1},~ \pi_{i2})$ represent the `unobserved' sensitivity and specificity respectively.  

Given study-specific sensitivity and specificity, two separate binomial distributions describe the distribution of true positives and true negatives among the diseased and the healthy individuals as follows
$$Y_{ij}~ |~ \pi_{ij}, ~\textbf{x}_i~ \sim~ bin(\pi_{ij},~ N_{ij}), i~=~1, ~\dots ~n, ~j~=~1, ~2$$

where $\textbf{x}_i$ generically denotes one or more covariates, possibly affecting $\pi_{ij}$. Equation ~\ref{eq:5} forms the higher level of the hierarchy and models the within-study variability. The second level of the hierarchy aims to model the between study variability of sensitivity and specificity while accounting for the inherent negative correlation thereof, with a bivariate distribution as follows
$$
\begin{pmatrix}
g(\pi_{i1})\\
g(\pi_{i2}) 
\end{pmatrix} \sim f(g(\pi_{i1}),~ g(\pi_{i2}))~ =~ f(g(\pi_{i1})) ~f(g(\pi_{i2})) ~c(F_1(g(\pi_{i1})),~ F_2(g(\pi_{i2}))),
$$		
where $g(.)$ denotes a transformation that might be used to modify the (0, 1) range to the whole real line. While it is critical to ensure that the studies included in the meta-analysis satisfy the specified entry criterion, there are study specific characteristics like different test thresholds and other unobserved differences that give rise to the second source of variability, the between-study variability. It is indeed the difference in the test thresholds between the studies that gives rise to the correlation between sensitivity and specificity. Including study level covariates allows us to model part of the between-study variability. The covariate information can and should be used to model the mean as well as the correlation between sensitivity and specificity.

In the next section we give more details on different bivariate distributions $f(g(\pi_{i1}),~g(\pi_{i2}))$ constructed using the logit or identity link function $g(.)$, different marginal densities and/or different copula densities $c$. We discuss their implications and demonstrate their application in meta-analysis of diagnostic accuracy studies.  An overview of suitable parametric families of copula for mixed models for diagnostic test accuracy studies was recently given by Nikoloupolous (2015). Here, we consider five copula functions which can model negative correlation.

#### Bivariate Gaussian copula
Given the density and the distribution function of the univariate and bivariate standard normal distribution with correlation parameter $\rho \in (-1, 1)$, the bivariate Gaussian copula function and density is expressed [@Meyer] as
$$C(u, ~v, ~\rho) = \Phi_2(\Phi^{-1}(u),~ \Phi^{-1}(v),~ \rho), $$
$$c(u, ~v, ~\rho) =~  \frac{1}{\sqrt{1~-~\rho^2}}  ~exp\bigg(\frac{2~\rho~\Phi^{-1}(u) ~\Phi^{-1}(v) - \rho^2~ (\Phi^{-1}(u)^2 + \Phi^{-1}(v)^2)}{2~(1 - \rho^2)}\bigg ) $$					

The logit transformation is often used in binary logistic regression to relate the probability of "success" (coded as 1, failure as 0) of the binary response variable with the linear predictor model that theoretically can take values over the whole real line. In diagnostic test accuracy studies, the `unobserved' sensitivities and specificities can range from 0 to 1 whereas their logits = $log ( \frac{\pi_{ij}}{1~-~ \pi_{ij}})$ can take any real value allowing to use the normal distribution as follows 
$$
logit (\pi_{ij}) ~\sim~ N(\mu_j, ~\sigma_j) ~<=> ~logit (\pi_{ij}) ~=~ \mu_j ~+~ \varepsilon_{ij},
$$
where, $\mu_j$ is a vector of the mean sensitivity and specificity for a study with zero random effects, and $\varepsilon_{i}$ is a vector of random effects associated with study $i$.
Now $u$ is the normal distribution function of $logit(\pi_{i1}$) with parameters $\mu_1$ and $\sigma_1$, \textit{v} is the normal distribution function of $logit(\pi_{i2})$ with parameters $\mu_2$ and $\sigma_2$,  $\Phi_2$ is the distribution function of a bivariate standard normal distribution with correlation parameter $\rho \in (-1, ~1)$ and $\Phi^{-1}$  is the quantile of the standard normal distribution. In terms of $\rho$, Kendall's tau is expressed as ($\frac{2}{\pi}$)arcsin$(\rho)$. 

With simple algebra the copula density with normal marginal distributions simplifies to 
$$
c(u, ~v, ~\rho) = \frac{1}{\sqrt{1 - \rho^2}} ~exp\bigg(\frac{1}{2 ~(1 ~-~ \rho^2)}~ \bigg( \frac{2~\rho~(x - \mu_1)~(y - \mu_2)}{\sigma_1 ~\sigma_2} ~-~ \rho^2 ~\bigg (\frac{{x ~-~ \mu_1}^2}{\sigma_1} ~+~ \frac{{y ~-~ \mu_2}^2}{\sigma_2}\bigg)\bigg ) \bigg ). $$

The product of the copula density above, the normal marginal of $logit(\pi_{i1}$) and $logit(\pi_{i2}$) form a bivariate normal distribution which characterize the model by [@Reitsma], [@Arends], [@Chu], and [@Rileyb], the so-called bivariate random-effects meta-analysis (BRMA) model, recommended as the appropriate method for meta-analysis of diagnostic accuracy studies. 
Study level covariate information explaining heterogeneity is introduced through the parameters of the marginal and the copula as follows 
$$\boldsymbol{\mu}_j = \textbf{X}_j\textbf{B}_j^{\top}. $$								
$\boldsymbol{X}_j$ is a $n \times p$ matrix containing the covariates values for the mean sensitivity($j = 1$) and specificity($j = 2$). For simplicity, assume that $\boldsymbol{X}_1 = \boldsymbol{X}_2 = \boldsymbol{X}$.   $\boldsymbol{B}_j^\top$ is a $p \times$ 1$ vector of regression parameters, and $p$ is the number of parameters.
By inverting the logit functions, we obtain
$$
\pi_{ij} = logit^{-1} (\mu_j + \varepsilon_{ij}).
$$									
Therefore, the meta-analytic sensitivity and specificity obtained by averaging over the random study effect, is given by, for $j = 1, 2$
$$
E(\pi_{j} )~=~E(logit^{-1} (\mu_j ~+~ \varepsilon_{ij}))
~=~\int_{-\infty}^{\infty}logit^{-1}(\mu_j ~+~ \varepsilon_{ij}) f(\varepsilon_{ij},~\sigma_j)~d\varepsilon_{ij},
$$								
assuming that $\sigma_1^2 > 0$ and $\sigma_2^2 > 0$. The integration above has no analytical expression and therefore needs to be numerically approximated and the standard are not easily available. Using MCMC simulation in the Bayesian framework the meta-analytic estimates can be easily computed as well as a standard error estimate and a credible intervals $E(\pi_{j})$ with minimum effort by generating predictions of the fitted bivariate normal distribution.

In the frequentist framework, it is more convenient however to use numerical averaging by sampling a large number $M$ of random-effects $\hat{\varepsilon}_{ij}$ from the fitted distribution and to estimate the meta-analytic sensitivity and specificity by, for $j = 1, 2$
$$
\hat{E}(\pi_{j}) = \frac{1}{M} \sum_{i~ =~ 1}^{M}logit^{-1}(\hat{\mu}_j + \hat{\varepsilon}_{ij}).
$$ 								
However, inference is not straightforward in the frequentist framework since the standard errors are not available. When $\varepsilon_{ij} = 0$, then
$$
E(\pi_{j} ~|~ \varepsilon_{ij} = 0) ~=~ logit^{-1}(\mu_j).
$$ 						
Inference for $E(\pi_{j} ~|~\varepsilon_{ij} ~=~ 0)$, as expressed in Equation~\ref{eq:14}, can be done in both Bayesian and frequentist framework.  The equation represents the mean sensitivity and specificity for a ``central" study with $\varepsilon_{ij} ~=~ 0$. Researchers often seem to confuse $E(\pi_{j} ~|~\varepsilon_{ij} ~=~ 0)$ with $E(\pi_{j})$ but due to the non-linear logit transformations, they are clearly not the same parameter. 

With the identity link function, no transformation on study-specific sensitivity and specificity is performed.  A natural choice for \textit{u} and \textit{v} would be beta distribution functions with parameters ($\alpha_1, ~\beta_1$) and ($\alpha_2, ~\beta_2$) respectively. Since $\pi_{ij} ~\sim ~ beta(\alpha_j, ~\beta_j)$, the meta-analytic sensitivity and specificity are analytically solved as follows
$$
E(\pi_{j} )
= \frac{\alpha_j}{\alpha_j + \beta_j},
$$

After reparameterising the beta distributions using the mean ($\mu_j ~=~ \frac{\alpha_j}{\alpha_j ~+~ \beta_j}$) and certainty ($\psi_j ~=~ \alpha_j ~+~ \beta_j$) or dispersion ($\varphi_j ~=~ \frac{1}{1~ +~ \alpha_j~+~\beta_j}$) parameters different link functions introduce covariate information to the mean, certainty/dispersion and association ($\rho$) parameters. A typical model parameterisation is 
$$
\boldsymbol{\mu}_j ~=~ logit^{-1}(\textbf{XB}_j^{\top}),  \\
\boldsymbol{\psi}_j ~=~ g(\textbf{WC}_j^{\top}),  \\
\boldsymbol{\alpha}_j ~=~ \boldsymbol{\mu}_{j}  ~\circ~ \boldsymbol{\psi}_{j}, \\
\boldsymbol{\beta}_j ~=~ ( \textbf{1} ~-~ \boldsymbol{\mu}_{j})  ~\circ~ \boldsymbol{\psi}_{j}, \\
\boldsymbol{\rho} ~=~ tanh(\textbf{ZD}_j^{\top}) ~=~ \frac{exp(2 \times \textbf{ZD}_j^{\top}) ~-~ 1}{expa(2\times \textbf{ZD}_j^{\top}) ~+~ 1}. $$								
$\boldsymbol{X, W}$ and $\boldsymbol{Z}$ are a $n \times p$ matrices containing the covariates values for the mean, dispersion and correlation which we will assume has similar information and denoted by $\boldsymbol{X}$ for simplicity purpose, $p$ is the number of parameters,  $\textbf{B}_j^{\top}$ , $\textbf{V}_j^{\top}$ and $\textbf{D}_j^{\top}$ are a $p \times 1$ vectors of regression parameters relating covariates to the mean, variance and correlation respectively. $g(.)$ is the log link to mapping $\textbf{XC}_j^{\top}$ to the positive real number line and $\circ$ is the Hadamard product. 

### Frank copula

This flexible copula in the so-called family of Archimedean copulas by [@Frank]. The functional form of the copula and the density is given by;
$$
C(F(\pi_{i1} ), ~F(\pi_{i2} ),\theta) ~=~ - \frac{1}{\theta} ~log \bigg [ 1 + \frac{(e^{-\theta ~F(\pi_{i1})} - 1)
    (e^{-\theta~ F(\pi_{i2})} - 1)}{e^{-\theta} - 1} \bigg ],	 \\				
c(F(\pi_{i1} ), ~F(\pi_{i2} ),\theta) ~=~  \frac {\theta~(1 - e^{-\theta})~e^{-\theta~(F(\pi_{i1}) ~ +~ F(\pi_{i2}))}}{[1 - e^{-\theta}-(1 - e^{-\theta ~F(\pi_{i1})})~(1 - e^{-\theta ~F(\pi_{i2})})]^2} . $$

Since $\theta \in \mathbb{R}$, both positive and negative correlation can be modelled, making this one of the more comprehensive copulas. When $\theta$ is 0, sensitivity and specificity are independent. For $\theta > 0$, sensitivity and specificity exhibit positive quadrant dependence and negative quadrant dependence when $\theta < 0$. The Spearman correlation $\rho_s$ and Kendall's tau $\tau_k$ can be expressed in terms of $\theta$ as
$$
\rho_s  = 1 - 12~\frac{D_2 (-\theta) ~-~ D_1(-\theta)}{\theta}, \\
\tau_k  = 1 + 4~ \frac{D_1(\theta) ~-~ 1}{\theta},
$$									
where $D_j(\delta)$ is the Debye function defined as
$$
D_j(\delta) ~=~ \frac{j}{\delta^j} \int_{\theta}^{\delta}\frac{t^j}{exp(t) ~-~ 1} ~dt, ~j~= 1, ~2.
$$

Covariate information is introduced in a similar manner as earlier. The identity link is used for the association parameter $\theta$. 

### Farlie-Gumbel-Morgenstern copula (FGM)

This popular copula by [@Farlie], [@Gumbel] and [@Morg] is defined as
$$
C(F(\pi_{i1} ), ~F(\pi_{i2}),~\theta) ~= ~F(\pi_{i1})~F(\pi_{i2})[1 ~+~ \theta~(1 - F(\pi_{i1}))~(1 ~-~ F(\pi_{i2}))], \nonumber \\				
c(F(\pi_{i1} ), ~F(\pi_{i2}),~\theta) ~=~ [1 ~+~ \theta~(2~F(\pi_{i1}) ~-~ 1)~(2~F(\pi_{i2}) ~-~ 1)].	$$				
Because $\theta \in (-1, 1)$, the Spearman correlation and Kendall's tau are expressed in terms of $\theta$ as $\theta/3$ and $2\theta/9$ respectively, making this copula only appropriate for data with weak dependence since $|\rho_s| ~\leq~ 1/3$. The logit link, log/identity link and Fisher's z transformation can be used to introduce covariate information in modelling the mean, dispersion and association parameter.

### Clayton copula

The Clayton copula function and density by [@Clayton] is defined as
$$
C(F(\pi_{i1}), ~F(\pi_{i2}),~\theta) = [F(\pi_{i1})^{- \theta} ~+~ F(\pi_{i2})^{-\theta} ~-~ 1]^{\frac{- 1}{\theta}}, \nonumber \\
c(F(\pi_{i1}), ~F(\pi_{i2}),~\theta) = (1 ~+~ \theta)~F(\pi_{i1})^{- (1 ~+~ \theta)} ~F(\pi_{i2})^{- (1 + \theta)}~[F(\pi_{i1} )^{- \theta}  + F(\pi_{i2})^{- \theta} ~-~ 1]^{\frac{- (2~\theta ~+~ 1)}{\theta}}. 	$$

Since $\theta \in (0, ~\infty)$, the Clayton copula typically models positive dependence; Kendall's tau equals $\theta/(\theta ~+~ 2)$. However, the copula function can be rotated by $90^\circ$ or $270^\circ$  to model negative dependence. The distribution and density functions following such rotations are given by
$$
C_{90} (F(\pi_{i1}), ~F(\pi_{i2}), ~\theta) ~=~ F(\pi_{i2}) ~-~ C(1 ~-~ F(\pi_{i1}), ~F(\pi_{i2}), ~\theta), \nonumber \\
c_{90} (F(\pi_{i1}),~ F(\pi_{i2}),~\theta) ~=~ (1 ~+~ \theta)(1 ~-~ F(\pi_{i1}))^{- (1 ~+~ \theta)}~F(\pi_{i2})^{- (1 ~+~ \theta)} ~[(1 - F(\pi_{i1}))^{-\theta} \nonumber \\ 
                                                                                                                                       +~ F(\pi_{i2})^{- \theta} - 1]^{\frac{- (2~\theta ~+~ 1)}{\theta}}, $$
and
$$
C_{270} (F(\pi_{i1}), ~F(\pi_{i2}), ~\theta) ~= ~F(\pi_{i1})- C(F(\pi_{i1}), ~1 ~-~ F(\pi_{i2}),\theta),  \nonumber \\ 					
c_{270} (F(\pi_{i1}), ~F(\pi_{i2}), ~\theta) ~=~(1 ~+~ \theta) ~F(\pi_{i1} )^{- (1 ~+~ \theta)}~ (1 ~-~ F(\pi_{i2} ))^{- (1 ~+~ \theta)}~[F(\pi_{i1})^{- \theta} \nonumber \\
                                                                                                                                           + (1 ~-~ F(\pi_{i2}))^{- \theta} ~-~ 1]^{\frac{- (2~\theta ~+~ 1)}{\theta}}. $$											
The logit, log/identity and log/identity links can be used to introduce covariate information in modelling the mean ($\mu_j$), certainty ($\psi_j$)/dispersion ($\varphi_j$) and association ($\theta$) parameters respectively.

Of course other copula functions that allow for negative association can be chosen. It is also an option to use known bivariate beta distributions. However, it is not always straightforward and analytically attractive to derive the corresponding copula function for all bivariate distributions. The use of existing bivariate beta distributions in meta-analysis of diagnostic accuracy studies has been limited because these densities model positive association( e.g., [@Libby], [@Olkina]), or both positive and negative association but over a restricted range( e.g., [@Sarmanov]). 

## Inference Framework and Software
Within the Bayesian framework, the analyst updates a prior opinion/information of a parameter based on the observed data whereas in the frequentist framework, the analyst investigates the behaviour of the parameter estimates in hypothetical repeated samples from a certain population. Due to its flexibility and use of MCMC simulations, complex modelling can often be implemented more easily within the Bayesian framework. By manipulating the prior distributions, Bayesian inference can circumvent identifiability problems whereas numerical approximation algorithms in frequentist inference without prior distributions can become stuck caused by identifiability problems. However, Bayesian methods typically require statistical expertise and patience because the MCMC simulations are computationally intensive. In contrast, most frequentist methods have been wrapped up in standard "procedures" that require less statistical knowledge and programming skills. Moreover frequentist methods are optimized with maximum likelihood estimation (MLE) that have much shorter run-times as opposed to MCMC simulations. `CopulaREMADA`[@copularemada] is such a MLE based `R` package.

### The CopulaDTA package
The `CopulaDTA` package is an extension of `rstan` [@rstan], the `R` interface to `Stan` [@stan] for diagnostic test accuracy data. Stan is a probabilistic programming language which has implemented Hamilton Monte Carlo(MHC) and uses No-U-Turn sampler (NUTS) [@Gelman]. The package facilitates easy application of complex models and their visualization within the Bayesian framework with much shorter run-times.

JAGS [@jags] is an alternative extensible general purpose sampling engine to Stan. Extending `JAGS` requires knowledge of `C++` to assemble a dynamic link library(DLL) module. From experience, configuring and building  the module is a daunting and tedious task especially in the Windows operation system. The above short-comings coupled with the fact that `Stan` tends to converge with fewer iterations even from bad initial values than `JAGS` made us prefer the `Stan` MCMC sampling engine.

The `CopulaDTA` package is available via the Comprehensive `R` Archive Network (CRAN) at https://CRAN.R-project.org/package=CopulaDTA. With a working internet connection, the `CopulaDTA` package is installed and loaded in `R` with the following commands


```r
#install.packages("CopulaDTA", dependencies = TRUE)
library(CopulaDTA)	
```

The `CopulaDTA` package provide functions to fit bivariate beta-binomial distributions constructed as a product of two beta marginal distributions and copula densities discussed earlier. The package also provides forest plots for a model with categorical covariates or with intercept only. Given the chosen copula function, a beta-binomial distribution is assembled up by the `cdtamodel` function which returns a `cdtamodel` object. The main function `fit` takes the `cdtamodel` object and fits the model to the given dataset and returns a `cdtafit` object for which `print`, `summary` and `plot` methods are provided for. 

### Model diagnostics
To assess model convergence, mixing and stationarity of the chains, it is necessary to check the potential scale reduction factor $\hat{R}$, effective sample size (ESS), MCMC error and trace plots of the parameters. When all the chains reach the target posterior distribution, the estimated posterior variance is expected to be close to the within chain variance such that the ratio of the two, $\hat{R}$ is close to 1 indicating that the chains are stable, properly mixed and likely to have reached the target distribution. A large $\hat{R}$ indicates poor mixing and that more iterations are needed. Effective sample size indicates how much information one actually has about a certain parameter. When the samples are auto correlated, less information from the posterior distribution of our parameters is expected than would be if the samples were independent. ESS close to the total post-warm-up iterations is an indication of less autocorrelation and good mixing of the chains. Simulations with higher ESS have lower standard errors and more stable estimates. Since the posterior distribution is simulated there is a chance that the approximation is off by some amount; the Monte Carlo (MCMC) error. MCMC error close to 0 indicates that one is likely to have reached the target distribution.

### Model comparison and selection
Watanabe-Alkaike Information Criterion (WAIC) [@Watanabe], a recent model comparison tool to measure the predictive accuracy of the fitted models in the Bayesian framework, will be used to compare the models. WAIC can be viewed as an improvement of Deviance Information Criterion(DIC) which, though popular, is known to be have some problems [@Plummer]. WAIC is a fully Bayesian tool, closely approximates the Bayesian cross-validation, is invariant to reparameterisation and can be used for simple as well as hierarchical and mixture models.

## Datasets
### Telomerase data}
[@Glas] systematically reviewed the sensitivity and specificity of cytology and other markers including telomerase for primary diagnosis of bladder cancer. They fitted a bivariate normal distribution to the logit transformed sensitivity and specificity values across the studies allowing for heterogeneity between the studies. From the included 10 studies, they reported that telomerase had a sensitivity and specificity of `0.75 [0.66, 0.74]` and `0.86 [0.71, 0.94]` respectively. They concluded that telomerase was not sensitive enough to be recommended for daily use. This dataset is available within the package and the following commands


```r
data(telomerase)
telomerase
```

loads the data into the R enviroment and generates the following output


```
##    ID TP  TN FN FP
## 1   1 25  25  8  1
## 2   2 17  11  4  3
## 3   3 88  31 16 16
## 4   4 16  80 10  3
## 5   5 40 137 17  1
## 6   6 38  24  9  6
## 7   7 23  12 19  0
## 8   8 27  18  6  2
## 9   9 14  29  3  3
## 10 10 37   7  7 22
```

`ID` is the study identifier, `DIS` is the number of diseased, `TP` is the number of true positives, `NonDis` is the number of healthy and `TN` is the number of true negatives. 

### ASCUS triage data
[@Arbyn] performed a Cochrane review on the accuracy of human papillomavirus testing and repeat cytology to triage of women with an equivocal Pap smear to diagnose cervical precancer. They fitted the BRMA model in `SAS` using `METADAS` on 10 studies where both tests were used. They reported absolute sensitivity of `0.909 [0.857, 0.944]` and `0.715 [0.629, 0.788]` for HC2 and repeat cytology respectively. The specificity was `0.607 [0.539, 0.68]` and `0.684 [0.599, 0.758]` for HC2 and repeat cytology respectively.  These data is used to demonstrate how the intercept only model is extended in a meta-regression setting. This dataset is also available within the package and the following commands


```r
data(ascus)
ascus
```

loads the data into the R enviroment and generates the following output



```
##    Test         StudyID  TP   FP  TN FN
## 1  RepC  Andersson 2005   6   14  28  4
## 2  RepC   Bergeron 2000   8   28  71  4
## 3  RepC Del Mistro 2010  20  191 483  7
## 4  RepC Kulasingam 2002  20   74 170  6
## 5  RepC     Lytwyn 2000   4   20  26  2
## 6  RepC      Manos 1999  48  324 570 15
## 7  RepC  Monsonego 2008  10   18 168 15
## 8  RepC      Morin 2001  14  126 214  5
## 9  RepC  Silverloo 2009  24   43 105 10
## 10 RepC    Solomon 2001 227 1132 914 40
## 11  HC2  Andersson 2005   6   17  25  4
## 12  HC2   Bergeron 2000  10   38  61  2
## 13  HC2 Del Mistro 2010  27  154 566  2
## 14  HC2 Kulasingam 2002  23  115 129  3
## 15  HC2     Lytwyn 2000   4   19  33  1
## 16  HC2      Manos 1999  58  326 582  7
## 17  HC2  Monsonego 2008  22  110  72  2
## 18  HC2      Morin 2001  17   88 253  2
## 19  HC2  Silverloo 2009  34   65  81  2
## 20  HC2    Solomon 2001 256 1050 984 11
```

`Test` is an explanatory variable showing the type of triage test, `StudyID` is the study identifier, `TP` is the number of true positives, `FP` is the number of false positives, `TN` is the number of true negatives, `FN` is the number of false negatives.













































