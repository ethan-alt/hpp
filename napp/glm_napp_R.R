library(cmdstanr)
model <- cmdstanr::cmdstan_model(
  'glm_napp.stan',
  include_paths = '~/Projects/Paper3/extensive_sims/napp'
)


get.dist.link = function(family) {
  fams  = c('binomial', 'poisson', 'gaussian', 'Gamma', 'inverse.gaussian')
  links = c(
    'identity', 'log', 'logit', 'inverse', 'probit', 'cauchit', 'cloglog'
    ,'sqrt', '1/mu^2'
  )
  fam.id  = which(family$family == fams)
  link.id = which(family$link == links)
  c(fam.id, link.id)
}



data.checks = function(
  formula, family, data, histdata, offset, offset0, check.hist = TRUE
) {
  if ( !(is.data.frame(data) ) )
    stop('data must be a data.frame')
  if ( any( is.na(data) ) )
    stop("data cannot contain missing values")
  if ( ! ( class(formula) == 'formula' ) )
    stop('formula must be of type "formula"')
  if ( class(family) != 'family')
    stop('family must be of type "family" (e.g., cannot be a character--use binomial() instead of "binomial"). See help(family)')
  if ( !formula.tools::is.two.sided(formula) )
    stop('formula must be two-sided')
  yname = formula.tools::lhs.vars(formula)
  if (length(yname) != 1)
    stop('formula must contain exactly 1 lhs variable name')
  varnames = all.vars(formula)
  if ( !( all( varnames %in% names(data) ) ) )
    stop("formula contains terms not in data")
  if (!is.null(offset))
    if(length(offset) != nrow(data))
      stop('offset and data must have equal lengths if offset is not NULL')
  if(check.hist) {
    if ( !(is.data.frame(histdata) ) )
      stop('histdata must be a data.frame')
    if ( any(is.na(histdata)) )
      stop("histdata cannot contain missing values")
    if ( !( all( varnames %in% names(histdata) ) ) )
      stop("formula contains terms not in histdata")
    if (!is.null(offset0))
      if(length(offset0) != nrow(histdata))
        stop('offset0 and histdata must have equal lengths if offset0 is not NULL')
  }
}



#'
#' Normalized asymptotic power prior (NAPP)
#'
#' Sample from the posterior distribution of a GLM using the normalized
#' asymptotic power prior. The regression coefficients and logarithm of the
#' dispersion parameter are a multivariate normal distribution with mean
#' equal to the maximum likelihood estimate of the historical data and
#' covariance matrix equal to \eqn{a0^{-1}} multiplied by the inverse Fisher
#' information matrix of the historical data, where a0 is the power prior
#' parameter (treated as random).
#'
#' @include data_checks.R
#'
#' @export
#'
#' @param formula      a two-sided formula giving the relationship between the response variable and covariates
#' @param family       an object of class `family`. See \code{\link[stats:family]{?stats::family}}
#' @param data         a `data.frame` giving the current data
#' @param histdata     a `data.frame` giving the historical data
#' @param a0.shape1    \code{shape1} parameter for a beta prior on the power prior parameter. When \code{a0.shape1 == 1} and \code{a0.shape2 == 1}, a uniform prior is used
#' @param a0.shape2    \code{shape1} parameter for a beta prior on the power prior parameter. When \code{a0.shape1 == 1} and \code{a0.shape2 == 1}, a uniform prior is used
#' @param offset       vector whose dimension is equal to the rows of the current data set giving an offset for the current data. Defaults to a vector of 0s
#' @param offset0      vector whose dimension is equal to the rows of the historical data set giving an offset for the historical data. Defaults to a vector of 0s
#' @param ...          arguments passed to [rstan::sampling()] (e.g. iter, chains).
#'
#' @return             an object of class `stanfit` giving posterior samples
#'
#' @examples
#' data(actg019)
#' data(actg036)
#' ## take subset for speed purposes
#' actg019 = actg019[1:100, ]
#' actg036 = actg036[1:50, ]
#' glm.napp(
#'   cd4 ~ treatment + age + race,
#'   family = poisson(), data = actg019, histdata = actg036,
#'   chains = 1, warmup = 500, iter = 1000
#' )
glm_napp = function(
  formula, family, data, histdata,
  a0.shape1    = 1.0,
  a0.shape2    = 1.0,
  offset       = NULL,
  offset0      = NULL,
  ...
) {
  y  = data[, all.vars(formula)[1]]
  y0 = histdata[, all.vars(formula)[1]]
  n  = length(y)
  n0 = length(y0)
  X  = model.matrix(formula, data)
  X0 = model.matrix(formula, histdata)
  fam.indx = get.dist.link(family)
  dist     = fam.indx[1]
  link     = fam.indx[2]

  ## Default offset is vector of 0s
  if ( is.null(offset) )
    offset = rep(0, n)
  if ( is.null(offset0) )
    offset0 = rep(0, n0)
  histdata$offset0 = offset0

  ## Fit enriched GLM to get hessian
  fit.glm =
    suppressWarnings(
      stats::glm(formula = formula, family = family, data = histdata, offset = offset0)
    )
  fit.glm    = enrichwith::enrich(fit.glm)
  theta.mean = fit.glm$coefficients
  theta.cov  = fit.glm$expected_information_mle

  if ( !( family$family %in% c('binomial', 'poisson') ) ) {
    ## theta = (beta, log(dispersion) )
    theta.mean = c(theta.mean, log(fit.glm$dispersion_mle))

    ## jacobian adjustment to fisher information for dispersion --> log dispersion
    theta.cov[nrow(theta.cov), nrow(theta.cov)] =
      fit.glm$dispersion_mle^2 * theta.cov[nrow(theta.cov), nrow(theta.cov)]
  }
  theta.cov = chol2inv(chol(theta.cov))

  ## perform data checks
  data.checks(
    formula, family, data, histdata, offset, offset0, check.hist = FALSE
  )


  standat = list(
    'n'            = n,
    'p'            = ncol(X),
    'ind_disp'     = ifelse(dist > 2, 1, 0),
    'y'            = y,
    'X'            = X,
    'theta_hat'    = theta.mean,
    'theta_cov'    = theta.cov,
    'a0_shape1'    = a0.shape1,
    'a0_shape2'    = a0.shape2,
    'dist'         = dist,
    'link'         = link,
    'offset'       = offset
  )

  ## fit model in cmdstanr --> convert to rstan
  fit = model$sample(data = standat, ...)
  # fit = rstan::read_stan_csv(fit$output_files())
  # 
  # ## rename parameters
  # newnames = names(fit.glm$coefficients)
  # if ( !( family$family %in% c('binomial', 'poisson') ) ) {
  #   newnames = c(newnames, 'dispersion')
  # }
  # fit@sim$fnames_oi[seq_along(newnames)] = newnames
  return(fit)
}









