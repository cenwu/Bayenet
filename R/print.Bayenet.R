#' print a Bayenet object
#'
#' Print a summary of a Bayenet object
#'
#' @param x Bayenet object.
#' @param digits significant digits in printout.
#' @param ... other print arguments.
#' @return No return value, called for side effects.
#' @usage \method{print}{Bayenet}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{Bayenet}}
#' @export
print.Bayenet=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nCoefficients:\n")
  print(x$coefficient, digits)
}



#' print a predict.Bayenet object
#'
#' Print a summary of a predict.Bayenet object
#'
#' @param x predict.Bayenet object.
#' @param digits significant digits in printout.
#' @param ... other print arguments.
#' @return No return value, called for side effects.
#' @usage \method{print}{Bayenet.pred}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{predict.Bayenet}}
#' @export
print.Bayenet.pred=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nPMSE:\n")
  print(x$error, digits)
  cat("\npredicted ", length(x$y.pred), " y (list component y.pred)", sep = "")
}

