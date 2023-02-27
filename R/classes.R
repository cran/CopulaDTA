#' @name cdtamodel-class
#' @title Class cdtamodel
#' @description A cdtamodel class in the CopulaDTA package.
#' @docType class
#' @slot copula copula function, 'fgm', 'gauss', 'c90', '270', or 'frank'.
#' @slot modelcode character with the model code as returned by the model function
#' @slot modelargs list containing control parameters for the prior distributions
#' @family cdta
#' @seealso \link{cdtamodel}
#' @export
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}

setClass(Class="cdtamodel",
          representation=representation(
              copula = 'character',
              modelcode = 'character',
              modelargs = "list"))

#' @name cdtafit-class
#' @title Class cdtafit
#' @description A cdtafit class in the CopulaDTA package.
#' @docType class
#' @slot data a data-frame with no missing values containing TP, TN, FP, FN, 'SID' and co-variables(if necessary).
#' @slot SID A string indicating the name of the column with the study identifier.
#' @slot copula copula function, 'fgm', 'gauss', 'c90', '270', or 'frank'.
#' @slot modelargs list containing control parameters for the prior distributions.
#' @slot fit an object of class stanfit returned by the function sampling.
#' @family cdta
#' @seealso \link{fit}
#' @export
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}
#' @importClassesFrom rstan stanmodel stanfit


setClass(Class="cdtafit",
         representation=representation(
            data='data.frame',
            SID = 'character',
            copula = 'character',
            modelargs = "list",
            fit='stanfit'))
