#' simulated data for demonstrating the features of Bayenet.
#'
#' Simulated gene expression data for demonstrating the features of Bayenet.
#'
#' @docType data
#' @keywords datasets
#' @name dat
#' @aliases dat  X Y clin coef
#' @usage data("dat")
#' @format dat consists of four components: X, Y, clin, coef. 
#' @details
#'
#' \strong{The data model for generating Y}
#' 
#' Use subscript \eqn{i} to denote the \eqn{i}th subject. Let \eqn{(Y_{i}, X_{i}, clin_{i})} (\eqn{i=1,\ldots,n}) be
#' independent and identically distributed random vectors. \eqn{Y_{i}} is a continuous response variable representing the
#' cancer outcome and disease phenotype. \eqn{X_{i}} is the \eqn{p}--dimensional vector of genetic factors. The clinical factors
#' is denoted as the \eqn{q}-dimensional vector \eqn{clin_{i}}.
#' The \eqn{\epsilon} follows some heavy-tailed distribution. Considering the following model:
#' \deqn{Y_{i} = \alpha_{0} + \sum_{k=1}^{q}\gamma_{k}C_{ik}+\sum_{j=1}^{p}\beta_{j}X_{ij}+\epsilon_{i},}
#' where \eqn{\alpha_{0}} is the intercept, \eqn{\gamma_{k}}'s and \eqn{\beta_{j}}'s are the regression coefficients corresponding to effects of clinical factors and genetic variants, respectively.
#' Denote \eqn{\gamma=(\gamma_{1}, \ldots, \gamma_{q})^{T}}, \eqn{\beta=(\beta_{1}, \ldots, \beta_{p})^{T}}. 
#' Then model can be written as 
#' \deqn{Y_{i} = C_{i}\gamma + X_{i}\beta + \epsilon_{i}.}
#'  
#' @examples
#' data(dat)
#' dim(X)
#' @seealso \code{\link{Bayenet}}
NULL