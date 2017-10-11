#' Trace plot using ggplot2.
#'
#' @param x An cdtafit object from \link{fit}.
#' @param ... additional options. See \link[rstan]{stan_trace} for more details.
#' @return A ggplot object of the parameters of the models mean structure.
#' @examples
#' \dontrun{
#' fit1 <- fit(model1,
#'                 SID='ID',
#'                 data=telomerase,
#'                 iter=2000,
#'                 warmup=1000,
#'                 thin=1,
#'                 seed=3)
#'
#' traceplot(fit1)
#'
#' traceplot(fit1) +
#' theme(axis.text.x = element_text(size=10, colour='black'),
#'       axis.text.y = element_text(size=10, colour='black'),
#'       axis.title.x = element_text(size=10, colour='black'),
#'       strip.text = element_text(size = 10, colour='black'),
#'       axis.title.y= element_text(size=10, angle=0, colour='black'),
#'       strip.text.y = element_text(size = 10, colour='black'),
#'       strip.text.x = element_text(size = 10, colour='black'),
#'       plot.background = element_rect(fill = "white", colour='white'),
#'       panel.grid.major = element_blank(),
#'       panel.background = element_blank(),
#'       strip.background = element_blank(),
#'       axis.line.x = element_line(color = 'black'),
#'       axis.line.y = element_line(color = 'black'))
#' }
#'@export

#' @author Victoria N Nyaga

traceplot.cdtafit <- function(x, ...){
    #open new window
    g <- rstan::stan_trace(x@fit, pars=c('MUse', 'MUsp'), ...)
    #draws <- as.array(x@fit, pars=c('MUse', 'MUsp'), ...)
    #g <- bayesplot::mcmc_trace(draws)
    return(g)
    if (grDevices::dev.interactive()) grDevices::dev.new()
    print(g)	}


