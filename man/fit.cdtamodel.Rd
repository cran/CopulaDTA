% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.R
\name{fit.cdtamodel}
\alias{fit.cdtamodel}
\title{Fit copula based bivariate beta-binomial distribution to diagnostic data.}
\usage{
fit.cdtamodel(
  cdtamodel,
  data,
  SID,
  cores = 3,
  chains = 3,
  iter = 6000,
  warmup = 1000,
  thin = 10,
  ...
)
}
\arguments{
\item{cdtamodel}{An object of cdtamodel class from \link{cdtamodel}.}

\item{data}{A data-frame with no missing values containing TP, TN, FP, FN, 'SID' and co-variables(if necessary).}

\item{SID}{A string indicating the name of the column with the study identifier.}

\item{cores}{A positive numeric values specifying the number of cores to use to execute parallel sampling. When the hardware has more at least 4 cores,
the default is 3 cores and otherwise 1 core.}

\item{chains}{A positive numeric value specifying the number of chains, default is 3.}

\item{iter}{A positive numeric value specifying the number of iterations per chain. The default is 6000.}

\item{warmup}{A positive numeric value (<iter) specifying the number of iterations to be discarded(burn-in/warm-up). The default is 1000.}

\item{thin}{A positive numeric value specifying the interval in which the samples are stored. The default is 10.}

\item{...}{Other optional parameters as specified in \link[rstan]{stan}.}
}
\value{
An object of cdtafit class.
}
\description{
Fit copula based bivariate beta-binomial distribution to diagnostic data.
}
\examples{
data(telomerase)
model1 <-  cdtamodel(copula = 'fgm')

model2 <- cdtamodel(copula = 'fgm',
               modelargs=list(param=2,
                              prior.lse='normal',
                              par.lse1=0,
                              par.lse2=5,
                              prior.lsp='normal',
                              par.lsp1=0,
                              par.lsp2=5))

model3 <-  cdtamodel(copula = 'fgm',
               modelargs = list(formula.se = StudyID ~ Test - 1))
\dontrun{
fit1 <- fit(model1,
                SID='ID',
                data=telomerase,
                iter=2000,
                warmup=1000,
                thin=1,
                seed=3)


fit2 <- fit(model2,
                SID='StudyID',
                data=ascus,
                iter=2000,
                warmup=1000,
                thin=1,
                seed=3)
}

}
\references{
{Nyaga VN, Arbyn M, Aerts M (2017). CopulaDTA: An R Package for Copula-Based Beta-Binomial Models for Diagnostic Test Accuracy
Studies in a Bayesian Framework. Journal of Statistical Software, 82(1), 1-27. doi:10.18637/jss.v082.c01}

{Agresti A (2002). Categorical Data Analysis. John Wiley & Sons, Inc.}

{Clayton DG (1978). A model for Association in Bivariate Life Tables and its Application in
Epidemiological Studies of Familial Tendency in Chronic Disease Incidence. Biometrika,65(1), 141-151.}

{Frank MJ (1979). On The Simultaneous Associativity of F(x, y) and x + y - F(x, y). Aequationes Mathematicae, pp. 194-226.}

{Farlie DGJ (1960). The Performance of Some Correlation Coefficients for a General Bivariate
Distribution. Biometrika, 47, 307-323.}

{Gumbel EJ (1960). Bivariate Exponential Distributions. Journal of the American Statistical Association, 55, 698-707.}

{Meyer C (2013). The Bivariate Normal Copula. Communications in Statistics - Theory and Methods, 42(13), 2402-2422.}

{Morgenstern D (1956). Einfache Beispiele Zweidimensionaler Verteilungen. Mitteilungsblatt furMathematische Statistik, 8, 23 - 235.}

{Sklar A (1959). Fonctions de Repartition a n Dimensions et Leurs Marges. Publications de l'Institut de Statistique de L'Universite de Paris, 8, 229-231.}
}
\author{
Victoria N Nyaga <victoria.nyaga@outlook.com>
}
