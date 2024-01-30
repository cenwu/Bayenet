#' @useDynLib Bayenet, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
#' fit a robust Bayesian elastic net variable selection model for genetic study.
#' @keywords models
#' @param X the matrix of predictors (genetic factors). Each row should be an observation vector. 
#' @param Y the continuous response variable. 
#' @param clin a matrix of clinical variables. Clinical variables are not subject to penalty. Clinical variables will be centered and a column of 1 will be added to the Clinical matrix as the intercept.
#' @param max.steps the number of MCMC iterations.
#' @param robust logical flag. If TRUE, robust methods will be used.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be used to shrink coefficients of irrelevant covariates to zero exactly.
#' @param penalty two choices are available. "lasso" for lasso penalty. "elastic net" for elastic net penalty.
#' @param debugging logical flag. If TRUE, progress will be output to the console and extra information will be returned.
#' @return an object of class `Bayenet' is returned, which is a list with component:
#' \item{posterior}{the posterior samples of coefficients from the MCMC.}
#' \item{coefficient}{the estimated value of coefficients.}
#' \item{burn.in}{the total number of burn-ins.}
#' \item{iterations}{the total number of iterations.}
#' \item{design}{the design matrix of all effects.}
#' 
#' @details Consider the data model described in "\code{\link{dat}}":
#' \deqn{Y_{i} = \alpha_{0} + \sum_{k=1}^{q}\gamma_{k}C_{ik}+\sum_{j=1}^{p}\beta_{j}X_{ij}+\epsilon_{i},}
#' where \eqn{\alpha_{0}} is the intercept, \eqn{\gamma_{k}}'s and \eqn{\beta_{j}}'s are the regression coefficients corresponding to effects of clinical factors and genetic variants, respectively.
#' 
#' When {penalty="elastic net"} (default), the elastic net penalty is adopted. If {penalty="lasso"}, the lasso penalty is used.
#' 
#' When {sparse=TRUE} (default), spike--and--slab priors are imposed to identify important main and interaction effects. If {sparse=FALSE}, Laplacian shrinkage will be used.
#'
#' When {robust=TRUE} (default), the distribution of \eqn{\epsilon_{i}} is defined as a Laplace distribution with density
#' \eqn{
#' f(\epsilon_{i}|\nu) = \frac{\nu}{2}\exp\left\{-\nu |\epsilon_{i}|\right\}
#' }, (\eqn{i=1,\dots,n}), which leads to a Bayesian formulation of LAD regression. If {robust=FALSE}, \eqn{\epsilon_{i}} follows a normal distribution.
#'
#' Both \eqn{X} and \eqn{clin} will be standardized before the generation of interaction terms to avoid the multicollinearity between main effects and interaction terms.
#'
#' Please check the references for more details about the prior distributions.
#' 
#' @references
#' Lu, X. and Wu, C. (2023). Bayesian quantile elastic net with spike-and-slab priors.
#' 
#' @examples
#' data(dat)
#'
#' max.steps=5000
#' fit= Bayenet(X, Y, clin, max.steps, penalty="lasso")
#' 
#' ## coefficients of parameters
#' fit$coefficient
#'
#' ## Estimated values of main G effects 
#' fit$coefficient$G
#' 
#' ## Estimated values of clincal effects 
#' fit$coefficient$clin
#' 
#' @export
Bayenet <- function(X, Y,clin, max.steps=10000, robust=TRUE, sparse=TRUE, penalty=c("lasso","elastic net"),debugging=FALSE)
{
  
  dat = DataMatrix(X, Y, clin, intercept=TRUE, debugging=FALSE)
  c=dat$c; g=dat$g; y=dat$y; beta_true=dat$coef
  n = dat$n; p= dat$p; q=ncol(c)
  
  G.names = dat$G.names
  clin.names = dat$clin.names
  
  if(robust){
    out = robust(X, Y, clin, max.steps, sparse, penalty,debugging=FALSE)
  }else{
    out = nonrobust(X, Y, clin, max.steps, sparse, penalty,debugging=FALSE)
  }
  
  BI=max.steps/2
  coeff.clin = apply(out$GS.alpha[-(1:BI),,drop=FALSE], 2, stats::median); 
  names(coeff.clin) = c(1,clin.names);
  
  coeff.G = apply(out$GS.beta[-(1:BI),,drop=FALSE], 2, stats::median); 
  names(coeff.G) = G.names;
  
  coefficient = list(clin=coeff.clin, G=coeff.G)
  out = list(GS.C=out$GS.alpha, GS.G=out$GS.beta)
  fit = list(posterior = out, coefficient=coefficient,burn.in = BI, iterations=max.steps, design=list(g,CLC=c))
  class(pred) = "Bayenet"
  fit
 
  
}
