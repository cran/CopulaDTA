#' Trace plot using ggplot2.
#'
#' @param object An object from \link{fitcopula}.
#' @param ... additional options. See \link[rstan]{stan_trace} for more details.
#' @return A ggplot trace plot of the parameters of the models mean structure.
#' @examples
#' \dontrun{
#' fit <- fitcopula(data=ascus,
#'          SID = "StudyID",
#'          formula.se= StudyID ~ Test,
#'          cores=3,
#'          seed=3,
#'          copula="fgm")
#'
#' tracecopula(fit)
#' }
#'@export
#' @author Victoria N Nyaga

tracecopula <- function(object, ...){
    #open new window
    g <- rstan::stan_trace(object$model, pars=c('MUse', 'MUsp'), ...)
    if (grDevices::dev.interactive()) grDevices::dev.new()
    print(g)
	}
