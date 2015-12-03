#' Fit copula based bivariate beta-binomial distribution to diagnostic data.
#'
#' @param data A data-frame with no missing values containg TP, TN, FP, FN, 'SID' and co-varaiables(if necessary).
#' @param SID A string indicating the name of the column with the study ID.
#' @param formula.se An optional object of class "formula": A symbolic description of a linear model to be fitted to mean E(x) of sensitivity in the logit scale.
#' the default (when no covariates are included) symbolic description is SID ~ 1 corresponds to the model formula E(x) = mu = exp(a)/(1 + exp(a)) where a is the intercept.
#' When the covariates are categorical and the relative measures are needed it is important to remove the interecept from the model to obtain meaningful parameters. EG for
#' a covariate 'Test' with two levels(A and B) and relative sensitivity of B versus A is needed, then the correct formula is SID ~ Test - 1 or SID ~ Test + 0. See \link[stats]{formula}.
#' For further information on interpretation of parameters in logistic regression see Agresti A(2002) Chapter 5.
#' @param formula.sp An optional object of class "formula": A symbolic description of a linear model to be fitted to specificity data.
#' By default the covariate information for sensitivity is used.
#' @param formula.omega An optional object of class "formula": A symbolic description of a linear model to be fitted to the copula function.
#' By default the covariate information for sensitivity is used.
#' @param copula a description of the copula function used to model the correlation between sensitivity and specificty.
#' This is a string naming the copula function. The choices are "fgm", "frank", "gauss", "c90" and "c270".
#' @param transform.omega A logical value indicating whether a constrained correlation parameter should be mapped into an non-constrained scale.
#' This applies to all the allowed copula functions except "frank". The default is TRUE.
#' @param param indication of the parameterisation used to map the marginal mean and precision/dispersion to the alpha and beta parameters of the beta distribution.
#' There are two choices: param=1 which uses \deqn{alpha = mu*phi, beta = (1 - mu)*phi} where \deqn{mu = alpha/(alpha + beta), 0<=mu<=1,} and
#' \deqn{phi = alpha + beta, phi >=0.}
#' param=2 uses \deqn{alpha = ((1 - phi)/phi)*mu}
#' \deqn{beta = ((1 - phi)/phi)*(1 - mu)} where \deqn{mu = alpha/(alpha + beta); 0<=mu<=1,} and
#' \deqn{phi = 1/( 1 + alpha + beta); 0<=phi<=1.}
#'@param prior.lse A description of prior distribution of the marginal mean sensitivity in the logit scale. The default is "normal" distribution.
#'For other distributions see stan documentation at \url{http://mc-stan.org/documentation/}.
#'@param par.lse1 A numeric value indicating the location of the prior distribution of the marginal mean sensitivity in the logit scale.
#'The default is 0 which implying a distribution centered around 0.5 in the 0-1 scale.
#'@param par.lse2 A numeric value indicating the spread(standard deviation) pf the prior distribution of the marginal mean sensitivity in the logit scale
#'and can be interpreted as the quantity of prior information.
#' vague and non-informative priors are specified by a distribution with large variance. The default is sd=10 implying that the variance is 100.
#'@param prior.lsp A description of prior distribution of the marginal mean specificity in the logit scale. The default is "normal" distribution.
#'@param par.lsp1 A numeric value indicating the location of the prior distribution of the marginal mean specificity in the logit scale.
#'The default is 0 which implying a distribution centered around 0.5 in the 0-1 scale.
#'@param par.lsp2 A numeric value indicating the spread(standard deviation) pf the prior distribution of the marginal mean specificity in the logit scale
#'and can be interpreted as the quantity of prior information.
#' vague and non-informative priors are specified by a distribution with large variance. The default is sd=10 implying that the variance is 100.
#'@param prior.omega A description of prior distribution of the correlation parameter(s). The default is "normal" distribution since "transform.omega=TRUE".
#'When "transform.omega=FALSE" the candidate prior distributions are U[-1, 1] for fgm and gaussian copulas, and half-cauchy(0, 2.5), gamma(0.001, 0.001) for the C90 and C270.
#'@param par.omega1 A numeric value indicating the location of the prior distribution of the correlation parameter(s). The default is 0.
#'@param par.omega2 A numeric value indicating the scale/spread(standard deviation) of the prior distribution of the correlation parameter(s). The default is sd=10.
#'@param chains A positive numeric value specifying the number of chains, default is 3.
#'@param iter A positive numeric value specifying the number of iterations per chain. The default is 6000.
#'@param warmup A positive numeric value (<iter) specifying the number of iterations to be discarded(burn-in/warm-up). The default is 1000.
#'@param thin A positive numeric value specifying the interval in which the samples are stored. The default is 10.
#'@param cores A positive numeric values specifying the number of cores to use to execute parallel sampling. When the hardware has more at least 4 cores,
#'the default is 3 cores and otherwise 1 core.
#'@param init One of digit 0, string "0" or "random", a function that returns a named list, or a list of named list.
#'"0": initialize all to be zero on the unconstrained support;
#'"random": randomly generated by Stan.
#'@param seed The random number generation seed for Stan. See \link[rstan]{stan} for more details.
#'@param ... Other oprtional parameters as specified in \link[rstan]{stan}.
#'@return An object which includes a S4 class stanfit. See \link[rstan]{stan} for more details.
#'
#'@examples
#' data(telomerase)
#' \dontrun{
#' fit1 <- fitcopula(data=telomerase,
#' SID = "ID",
#' copula="fgm",
#' iter = 4000,
#' warmup = 1000,
#' thin = 10,
#' seed=3,
#' cores=1)
#
#' data(ascus)
#'
#' fit2 <- fitcopula(data=ascus,
#' SID = "StudyID",
#' copula="fgm",
#' formula.se= StudyID ~ Test - 1,
#' seed=3)
#' }
#'
#'@references {Agresti A (2002). Categorical Data Analysis. John Wiley & Sons, Inc.}
#'@references {Clayton DG (1978). A model for Association in Bivariate Life Tables and its Application in
#'Epidemiological Studies of Familial Tendency in Chronic Disease Incidence. Biometrika,65(1), 141-151.}
#'@references {Frank MJ (1979). On The Simultaneous Associativity of F(x, y) and x + y - F(x, y). Aequationes Mathematicae, pp. 194-226.}
#'@references {Farlie DGJ (1960). The Performance of Some Correlation Coefficients for a General Bivariate
#'Distribution. Biometrika, 47, 307-323.}
#'@references {Gumbel EJ (1960). Bivariate Exponential Distributions. Journal of the American Statistical Association, 55, 698-707.}
#'@references {Meyer C (2013). The Bivariate Normal Copula. Communications in Statistics - Theory and Methods, 42(13), 2402-2422.}
#'@references {Morgenstern D (1956). Einfache Beispiele Zweidimensionaler Verteilungen. Mitteilungsblatt furMathematische Statistik, 8, 23 - 235.}
#'@references {Sklar A (1959). Fonctions de Repartition a n Dimensions et Leurs Marges. Publications de l'Institut de Statistique de L'Universite de Paris, 8, 229-231.}
#'@export
#'@importFrom rstan sampling
#'@importFrom rstan stan_model
#' @author Victoria N Nyaga
fitcopula <- function(data,
                      SID,
					  copula,
					  formula.se=NULL,
					  formula.sp=NULL,
					  formula.omega=NULL,
					  transform.omega=TRUE,
					  param=1,
					  prior.lse='normal',
					  par.lse1=0,
					  par.lse2=10,
					  prior.lsp='normal',
					  par.lsp1=0,
					  par.lsp2=10,
					  prior.omega=NULL,
					  par.omega1=NULL,
					  par.omega2=NULL,
					  chains=NULL,
					  iter = 6000,
					  warmup = 1000,
					  thin = 10,
					  cores=NULL,
					  init = "random",
					  seed = sample.int(.Machine$integer.max, 1), ...) {

#==================================Build model piece by piece =======================
pt1 <- "functions{"

#=================================Choice model ======================================#
GAUSS <- "\n\treal PhiInv(real u){\n\t\t
		real q;
		real v;
		real r;
		real val;

		q <- u - 0.5;
		if (0.075 <= u <= 0.925){
			v <- u - 0.5;
			r <- .180625 - v*v;
			val <-    q * (((((((r * 2509.0809287301226727 +
							   33430.575583588128105) * r + 67265.770927008700853) * r +
							 45921.953931549871457) * r + 13731.693765509461125) * r +
						   1971.5909503065514427) * r + 133.14166789178437745) * r +
						 3.387132872796366608)/
				(((((((r * 5226.495278852854561 +
						 28729.085735721942674) * r + 39307.89580009271061) * r +
					   21213.794301586595867) * r + 5394.1960214247511077) * r +
					 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
		}
		else {
			if (u < 0.5){
				r <- u;
			}
			else {
				r <- 1 - u;
			}
			r <- sqrt(-log(r));
			if (r<= 5) {
				r <- r - 1.6;
				val <- (((((((r * 7.7454501427834140764e-4 +
							.0227238449892691845833) * r +  .24178072517745061177) *r +
							1.27045825245236838258) * r + 3.64784832476320460504) * r +
						5.7694972214606914055) *r +      4.6303378461565452959) * r +
						   1.42343711074968357734)/
					(((((((r *1.05075007164441684324e-9 + 5.475938084995344946e-4) *r +
						  .0151986665636164571966) * r + .14810397642748007459) * r +
						  .68976733498510000455) * r + 1.6763848301838038494) * r +
						2.05319162663775882187) * r + 1);
			}
			else {
				r <- r - 5;
				val <- (((((((r * 2.01033439929228813265e-7    + 2.71155556874348757815e-5) * r +
								.0012426609473880784386) * r + .026532189526576123093) * r +
								.29656057182850489123) * r + 1.7848265399172913358) * r +
							  5.4637849111641143699) *  r + 6.6579046435011037772)/ (((((((
							r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7)*r +
							1.8463183175100546818e-5) * r + 7.868691311456132591e-4) * r +
						.0148753612908506148525) * r + .13692988092273580531) * r +
						.59983220655588793769) * r + 1);
			}
			if(q < 0) val <- -val;
		 }
		return(val);
	}
	real gaussian_log(matrix p_i, matrix alpha, matrix beta, vector omega){
		real f1;
		real f2;
		vector[rows(p_i)] f3;
		vector[rows(p_i)] f4;
		int r;

		r <- rows(p_i);
		f1 <- beta_log(col(p_i, 1), col(alpha,1), col(beta,1));
		f2 <- beta_log(col(p_i, 2), col(alpha,2), col(beta,2));
		for (i in 1:r){
			f3[i] <- (1/sqrt(1 - omega[i]^2))*exp((2*omega[i]*PhiInv(beta_cdf(p_i[i, 1], alpha[i,1], beta[i,1]))*PhiInv(beta_cdf(p_i[i, 2], alpha[i,2], beta[i,2])) -
			omega[i]^2*(PhiInv(beta_cdf(p_i[i, 1], alpha[i,1], beta[i,1]))^2 + PhiInv(beta_cdf(p_i[i, 2], alpha[i,2], beta[i,2]))^2))/
			(2*(1 - omega[i]^2)));
		}
		return (f1 + f2 + sum(log(f3)));
	}
"

FRANK <- "\n\treal frank_log(matrix p_i, matrix alpha, matrix beta, vector omega){
	real f1;
	real f2;
	int r;
	vector[rows(p_i)] f3;
	vector[rows(p_i)] f4;
	r <- rows(p_i);

	f1 <- beta_log(col(p_i, 1), col(alpha,1), col(beta, 1));
	f2 <- beta_log(col(p_i, 2), col(alpha, 2), col(beta, 2));
	for (i in 1:r){
		f3[i] <- (omega[i]*(1 - exp(-omega[i]))*exp(-omega[i]*(beta_cdf(p_i[i,1], alpha[i,1], beta[i,1]) + beta_cdf(p_i[i,2], alpha[i,2], beta[i,2]))));
		f4[i] <- ((1 - exp(-omega[i])) - (1 - exp(-omega[i]*beta_cdf(p_i[i,1], alpha[i,1], beta[i,1])))*(1 - exp(-omega[i]*beta_cdf(p_i[i,2], alpha[i,2], beta[i,2]))));
	}
	return (f1 + f2 + sum(log((f3)./((f4).*(f4)))));
}
"
FGM <- "\n\treal fgm_log(matrix p_i, matrix alpha, matrix beta, vector omega){
	real f1;
	real f2;
	vector[rows(p_i)] f3;
    int r;

    r <- rows(p_i);

	f1 <- beta_log(col(p_i, 1), col(alpha, 1), col(beta, 1));
	f2 <- beta_log(col(p_i, 2), col(alpha, 2), col(beta, 2));
	for (i in 1:r){
		f3[i] <- log(1 + omega[i]*(1 - 2*beta_cdf(p_i[i,1], alpha[i,1], beta[i,1]))*(1 - 2*beta_cdf(p_i[i,2], alpha[i,2], beta[i,2])));
	}
	return (f1 + f2 + sum(f3));
	}
"
C90 <- "\n\treal clayton90_log(matrix p_i, matrix alpha, matrix beta, vector omega){
	real f1;
	real f2;
	vector[rows(p_i)] f3;
	vector[rows(p_i)] f4;
	int r;
	vector[rows(p_i)] powr;

	r <- rows(p_i);
	f1 <- beta_log(col(p_i, 1), col(alpha,1), col(beta,1));
	f2 <- beta_log(col(p_i, 2), col(alpha,2), col(beta,2));
	powr <- -(2*omega + 1)./(omega);

	for (i in 1:r){
		f3[i] <- (1 + omega[i])*((1 - beta_cdf(p_i[i,1], alpha[i,1], beta[i,1]))^(-(1 + omega[i])))*(beta_cdf(p_i[i,2], alpha[i,2], beta[i,2])^(-(1 + omega[i])))*
		(((1 - beta_cdf(p_i[i,1], alpha[i,1], beta[i,1]))^(-omega[i]) + beta_cdf(p_i[i,2], alpha[i,2], beta[i,2])^(-omega[i]) - 1)^powr[i]);
	}
	return (f1 + f2 + sum(log(f3)));
	}
"
C270 <- "\n\treal clayton270_log(matrix p_i, matrix alpha, matrix beta, vector omega){
	real f1;
	real f2;
	vector[rows(p_i)] f3;
	vector[rows(p_i)] f4;
	int r;
	vector[rows(p_i)] powr;

	r <- rows(p_i);
	f1 <- beta_log(col(p_i, 1), col(alpha,1), col(beta,1));
	f2 <- beta_log(col(p_i, 2), col(alpha,2), col(beta,2));
	powr <- -(2*omega + 1)./(omega);

	for (i in 1:r){
		f3[i] <- (1 + omega[i])*(beta_cdf(p_i[i,1], alpha[i,1], beta[i,1])^(-(1 + omega[i])))*((1 - beta_cdf(p_i[i,2], alpha[i,2], beta[i,2]))^(-(1 + omega[i])))*
		((beta_cdf(p_i[i,1], alpha[i,1], beta[i,1])^(-omega[i]) + (1 - beta_cdf(p_i[i,2], alpha[i,2], beta[i,2]))^(-omega[i]) - 1)^powr[i]);
	}
	return (f1 + f2 + sum(log(f3)));
	}
"

if (copula=="gauss"){
		copkies<-GAUSS
	} else if (copula=="frank") {
		 copkies<-FRANK
	} else if (copula=="fgm") {
		 copkies<-FGM
	} else if (copula=="c90") {
		 copkies<-C90
	} else if (copula=="c270") {
		 copkies<-C270
	} else {
		stop("Invalid copula chosen")}

#========================== Data ============================================#

df <- prep.data(data,
                SID = SID,
                formula.se=formula.se,
                formula.sp=formula.sp,
                formula.omega=formula.omega)

datalist <- list(
	Ns = df$Ns,
	Npse = df$Npse,
	Npsp = df$Npsp,
	Npomega = df$Npomega,
	tp = df$data$TP,
	dis = df$data$TP + df$data$FN,
	tn = df$data$TN,
	nondis = df$data$TN + df$data$FP,
	xse = df$XSE,
	xsp = df$XSP,
	xomega = df$XOMEGA)

dat <- "\n}\n data{
		int<lower=0> Ns;
		int<lower=0> tp[Ns];
		int<lower=0> dis[Ns];
		int<lower=0> tn[Ns];
		int<lower=0> nondis[Ns];

		int<lower=0> Npse;
		matrix<lower=0>[Ns,Npse] xse;

		int<lower=0> Npsp;
		matrix<lower=0>[Ns,Npsp] xsp;

		int<lower=0> Npomega;
		matrix<lower=0>[Ns,Npomega] xomega;\n}"

#======================= Parameters ====================================#
params <- "\n parameters{
		vector[Npse] betamuse;
		vector[Npsp] betamusp;
		vector[Npse] betaphise;
		vector[Npsp] betaphisp;
		matrix<lower=0, upper=1>[Ns,2] p_i;
		"

if(copula=="frank" | transform.omega==TRUE){
	betaomega <- "vector[Npomega] betaomega;"
} else if (copula=="gauss" |copula=="fgm" ){
	 betaomega <- "vector<lower=-1, upper=1>[Npomega] betaomega;"

} else if (copula=="c90" |copula=="c270" ){
	 betaomega <- "vector<lower=0>[Npomega] betaomega;"}

#======================= Transformed Parameters ====================================#
transf_params.pt1 <- "\n } \n transformed parameters{
		matrix<lower=0, upper=1>[Ns,2] mui;
		vector<lower=0, upper=1>[Ns] musei;
		vector<lower=0, upper=1>[Ns] muspi;
		vector[Npse] MUse;
		vector[Npsp] MUsp;
		vector[Npse] RRse;
		vector[Npsp] RRsp;
		matrix<lower=0>[Ns,2] alpha;
		matrix<lower=0>[Ns,2] beta;"

transf_params.pt2 <- "\n\n\t\tmusei <- exp(xse*betamuse)./(1 + exp(xse*betamuse));
		muspi <- exp(xsp*betamusp)./(1 + exp(xsp*betamusp));
		mui <- append_col(musei, muspi);

		MUse <- exp(betamuse)./(1 + exp(betamuse));
		MUsp <- exp(betamusp)./(1 + exp(betamusp));

		RRse[1] <- 1;
		for (i in 2:Npse)
			RRse[i] <- MUse[i]/MUse[1];

		RRsp[1] <- 1;
		for (i in 2:Npsp)
			RRsp[i] <- MUsp[i]/MUsp[1];"
if(param==1){
phi.pt1 <- "\n\t\tvector<lower=0>[Ns] phisei; \n\t\tvector<lower=0>[Ns] phispi; \n\t\tmatrix<lower=0>[Ns,2] phi;\n"
phi.pt2 <- "\n\t\tphisei <- exp(xse*betaphise);
		phispi <- exp(xsp*betaphisp);
		phi <- append_col(phisei, phispi);
		alpha <- (mui).*phi;
		beta <- (1 - mui).*phi;\n"
}else{

phi.pt1 <- "\n\t\tvector<lower=0, upper=1>[Ns] phisei; \n\t\tvector<lower=0, upper=1>[Ns] phispi; \n\t\tmatrix<lower=0, upper=1>[Ns,2] phi;\n"
phi.pt2 <- "\n\t\tphisei <- exp(xse*betaphise)./(1 + exp(xse*betaphise));
		phispi <- exp(xsp*betaphisp)./(1 + exp(xsp*betaphisp));
		phi <- append_col(phisei, phispi);
		alpha <- ((1 - phi)./(phi)).*(mui);
		beta <- ((1 - phi)./(phi)).*(1 - mui);\n"
}

if (copula=="frank"){
	omega.pt1 <- "\t\t vector[Ns] omega;"
	omega.pt2 <- "\n\t\t omega <- xomega*betaomega;"
} else if (copula=="gauss"|copula=="fgm"){
	   if(transform.omega==TRUE){
			omega.pt1 <- "\t\tvector[Ns] omegat;\n\t\tvector<lower=-1, upper=1>[Ns] omega;\n\t\tvector[Npomega] betaomegat;"
			omega.pt2 <- "\n\t\tomegat <- xomega*betaomega; \n\t\tfor (s in 1:Ns) \n\t\t\tomega[s] <- tanh(omegat[s]);\n\t\tfor (o in 1:Npomega) \n\t\t\tbetaomegat[o] <- tanh(betaomega[o]);"
		}else{
			omega.pt1 <- "\t\tvector<lower=-1, upper=1>[Ns] omega;"
			omega.pt2 <- "\n\t\tomega <- xomega*betaomega;"
	}
} else {
	if(transform.omega==TRUE){
		omega.pt1 <- "\t\tvector[Ns] omegat;\n\t\tvector<lower=0>[Ns] omega;\n\t\tvector[Npomega] betaomegat;"
		omega.pt2 <-"\n\t\tomegat <- xomega*betaomega; \n\t\tfor (s in 1:Ns) \n\t\t\tomega[s] <- exp(omegat[s]);\n\t\tfor (o in 1:Npomega) \n\t\t\tbetaomegat[o] <- exp(betaomega[o]);"
	}else{
		omega.pt1 <- "\t\t vector<lower=0>[Ns]omega;"
		omega.pt2 <- "\n\t\t omega <- xomega*betaomega;"
	}
}

# if (copula=="frank"){
# 	ktau.pt1 <- "\n"
# } else{
	ktau.pt1 <- "\n\t\tvector[Npomega] ktau;"
#}

if (copula=="frank"){
    ktau.pt2 <- "\n\t\tfor (o in 1:Npomega) \n\t\t\tktau[o] <- 1;"
} else if (copula=="gauss"){
	if(transform.omega==TRUE){
		ktau.pt2 <- "\n\t\tfor (o in 1:Npomega) \n\t\t\tktau[o] <- (2/pi())*asin(betaomegat[o]);"
	}else{
		ktau.pt2 <- "\n\t\tfor (o in 1:Npomega) \n\t\t\tktau[o] <- (2/pi())*asin(betaomega[o]);"
	}
} else if(copula=="fgm"){
   if(transform.omega==TRUE){
		ktau.pt2 <- "\n\t\tfor (o in 1:Npomega) \n\t\t\tktau[o] <- (2*betaomegat[o])/9;"
	}else{
		ktau.pt2 <- "\n\t\tfor (o in 1:Npomega) \n\t\t\tktau[o] <- (2*betaomega[o])/9;"
	}
} else {
	if(transform.omega==TRUE){
		ktau.pt2 <- "\n\t\tfor (o in 1:Npomega) \n\t\t\tktau[o] <- -(betaomegat[o])/(betaomegat[o] + 2);"
	}else{
		ktau.pt2 <- "\n\t\tfor (o in 1:Npomega) \n\t\t\tktau[o] <- -(betaomega[o])/(betaomega[o] + 2);"
	}
}

#======================================Priors ===========================================#
pt2 <- "\n}\n model{\n\t #priors \n"

priorse <- paste('\t betamuse ~ ', prior.lse,'(', par.lse1, ', ', par.lse2, ');\n', sep='')
priorsp <- paste('\t betamusp ~ ', prior.lsp,'(', par.lsp1, ', ', par.lsp2, ');\n', sep='')

if(copula=="frank"){if(is.null(prior.omega)){prioromega<-paste('\t betaomega ~ ', "normal",'(', 0, ', ', 10, ');\n', sep='')} else{
	prioromega<-paste('\t betaomega ~ ', prior.omega,'(', par.omega1, ', ', par.omega2, ');\n', sep='')

}} else {
	if(transform.omega==FALSE){if(is.null(prior.omega)){prioromega<-paste('\t betaomega ~ ', "normal",'(', 0, ', ', 10, ');\n', sep='')} else{
		prioromega<-paste('\t betaomega ~ ', prior.omega,'(', par.omega1, ', ', par.omega2, ');\n', sep='')

}} else {
	if(is.null(prior.omega)){prioromega<-paste('\t betaomega ~ ', "normal",'(', 0, ', ', 10, ');\n', sep='')} else {
		prioromega=paste('\t betaomega ~ ', prior.omega,'(', par.omega1, ', ', par.omega2, ');\n', sep='')}}}


if (copula=="gauss"){
	priorp <- "\n\t p_i ~ gaussian(alpha, beta, omega);\n"
} else if (copula=="frank"){
	priorp <- "\n\t p_i ~ frank(alpha, beta, omega);\n"
} else if (copula=="fgm"){
	priorp <- "\n\t p_i ~ fgm(alpha, beta, omega);\n"
} else if (copula=="c90"){
	priorp <- "\n\t p_i ~ clayton90(alpha, beta, omega);\n"
} else if(copula=="c270"){
	priorp <- "\n\t p_i ~ clayton270(alpha, beta, omega);\n"
}

#===============================Likelihood=============================================#

lik <- "\n\t tp ~ binomial(dis,col(p_i,1)); \n\t tn ~ binomial(nondis, col(p_i, 2)); \n }"

#===================================Generated Quantities=============================#
GQ <- "\n generated quantities{
	vector[Ns*2] loglik;
	for (i in 1:Ns)
		loglik[i] <- binomial_log(tp[i], dis[i], p_i[i,1]);

	for (i in (Ns+1):(2*Ns))
		loglik[i] <- binomial_log(tn[i-Ns], nondis[i-Ns], p_i[i-Ns,2]);\n}"

model <- paste(pt1,
		  copkies,
		  dat,
		  params,
		  betaomega,
		  transf_params.pt1,
		  phi.pt1,
		  omega.pt1,
		  ktau.pt1,
		  transf_params.pt2,
		  phi.pt2,
		  omega.pt2,
		  ktau.pt2,
		  pt2,
		  priorse,
		  priorsp,
		  prioromega,
		  priorp,
		  lik,
		  GQ,
		  sep='')
#=================================Run the model=========================================#
#Check available cores, if more than 4, use
if((parallel::detectCores() > 3) & is.null(cores)) {
	cores <- 3
} else if((parallel::detectCores() < 3) & is.null(cores)) {
	cores < -1
}

if(is.null(chains)) chains <- 3

sm <- rstan::stan_model(model_code = model)

mod <- rstan::sampling(object = sm,
             data=datalist,
             chains = chains,
             iter = iter,
             warmup = warmup,
             thin = thin,
             seed = seed,
             cores=cores, ...)
#=======================================FINISH=======================================#
    return(list(model=mod, copula=copula, data=data, SID=SID, chains=chains, iter=iter,
                warmup=warmup, thin=thin, formula.se=formula.se, formula.sp=formula.sp,
                formula.omega=formula.omega))
}
