#' @export
## ---- summary.ame----
summary.ame <- function(object, ...) {

  # summarize results for continuous variables
  if (!is.null(object$continuous)) {

    ame_mean <- sapply(object$continuous, mean)
    ame_sd   <- sapply(object$continuous, sd)
    ame_zval <- ame_mean/ame_sd
    ame_pval <- 1.96 * pnorm(-abs(ame_zval))

    res_continuous <- cbind(ame_mean, ame_sd, ame_zval, ame_pval)

  } else {

    res_continuous <- NULL

  }

  # summarize results for factor variables
  if (!is.null(object$factor)) {

    fd_mean <- sapply(object$factor, mean)
    fd_sd   <- sapply(object$factor, sd)
    fd_zval <- fd_mean/fd_sd
    fd_pval <- 1.96 * pnorm(-abs(fd_zval))

    res_factor <- cbind(fd_mean, fd_sd, fd_zval, fd_pval)

  } else {

    res_factor <- NULL

  }

  # put results into list
  res_list <- list(continuous = res_continuous, factor = res_factor)
  attr(res_list, "continuous_variables") <- attr(object$continuous, "variable_names")
  attr(res_list, "factor_variables") <- attr(object$factor, "variable_names")
  attr(res_list, "n_sim") <- attr(object, "n_sim")

  class(res_list) <- "summary.ame"

  return(res_list)
}

#' @export
# ---- print.summary.ame ----
print.summary.ame <- function(x, ...) {

  if (is.null(x$continuous)) {
    res <- data.frame(x$factor)
  } else if (is.null(x$factor)) {
    res <- data.frame(x$continuous)
  } else {
    res <- data.frame(do.call("rbind", x))
  }

  names(res) <- c("dydx", "Std. Error", "z value", "Pr(>|z|)")

  c_vars <- paste(attr(x, "continuous_variables"), sep = " ", collapse = ", ")
  f_vars <- paste(attr(x, "factor_variables"), sep = " ", collapse =", ")

  # print results to screen
  cat("Continuous variables:", c_vars, "\n", sep = " ")
  cat("Factor variables:    ", f_vars, "\n\n", sep = " ")
  printCoefmat(res, signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE, signif.legend = TRUE, ...)
  cat("\n")
  cat("Note: For a factor variable f, dydx corresponds to the first difference E(Y|f_i) - E(Y|f_0)")

}
