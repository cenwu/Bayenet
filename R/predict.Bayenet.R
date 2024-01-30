#' make predictions from a Bayenet object
#'
#' make predictions from a Bayenet object
#' 
#' @param object Bayenet object.
#' @param X.new a matrix of new values for X at which predictions are to be made.
#' @param clin.new a vector or matrix of new values for clin at which predictions are to be made.
#' @param Y.new a vector of the response of new observations. If provided, the prediction error will be computed based on Y.new.
#' @param ... other predict arguments
#' @details X.new must have the same number of columns as X used for fitting the model. If clin was provided when fit the model, clin.new
#' must not be NULL, and vice versa. The predictions are made based on the posterior estimates of coefficients in the Bayenet object.
#' Note that the effects of clinical factors are not subject to selection.
#'
#' If Y.new is provided, the prediction error will be computed. For robust methods, the prediction mean absolute deviations (PMAD) will be computed.
#' For non-robust methods, the prediction mean squared error (PMSE) will be computed.
#'
#' @return  an object of class `Bayenet.pred' is returned, which is a list with components:
#' \item{error}{prediction error. error is NULL is Y.new=NULL.}
#' \item{y.pred}{predicted values of the new observations.}
#'
#' @rdname predict.Bayenet
#' @seealso \code{\link{Bayenet}}
#'
#@examples
#data(dat)
#test=sample((1:nrow(X)), floor(nrow(X)/5))
#fit=Bayenet(X[-test,], Y[-test], clin[-test,], max.steps=5000,penalty="lasso")
#predict.Bayenet(fit, X[test,], clin[test,], Y[test,])
#'
#' @export
predict.Bayenet=function(object, X.new, clin.new, Y.new,...){

  intercept = TRUE
  dat = DataMatrix(X.new, Y.new, clin.new, intercept)
  c=dat$c; g=dat$g; y=dat$y; beta_true=dat$coef
  n = dat$n; p= dat$p; q=ncol(c)

  coeff = object$coefficient$G
  coeff.clc = object$coefficient$clin

  y.pred = g %*% coeff + c %*% coeff.clc
  error = NULL

  if(inherits(object, "RBVS")){
    error = sum(abs(Y.new - y.pred))/length(Y.new)
    # error.type = "PMAD"
    names(error) = "PMAD"
  }else{
    error = sum((Y.new - y.pred)^2)/length(Y.new)
    # error.type = "PMSE"
    names(error) = "PMSE"
  }

  pred = list(error=error, y.pred=y.pred)
  class(pred) = "Bayenet.pred"
  pred
}


