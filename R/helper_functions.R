# helper-functions
## ---- coefficient/covariate names ----
get_term_names <- function(coefnames) {
  name <- gsub("(", "", coefnames, fixed = TRUE)
  name <- gsub(")", "", name, fixed = TRUE)

  # get coefficient and covariate names
  betas <- paste("beta_", gsub(":", "X", name), sep = "")
  covariates <- gsub(":", " * ", name)

  term_names <- list(betas = betas, covariates = covariates)
  return(term_names)
}

## ---- obain term names in partial derivatives ----
get_deriv_terms <- function(d, dydx_vars) {
  deriv_terms <- list()
  for (i in 1:length(dydx_vars)) {
    d_n <- gsub(" ", "", as.character(d[[i]]))
    d_n <- strsplit(d_n, split = c("+"), fixed = TRUE)
    deriv_terms[[i]] <- lapply(d_n, function(x) strsplit(x, split = c("*"), fixed = TRUE))
  }
  deriv_terms <- unique(unlist(deriv_terms))
  return(deriv_terms)
}

## ---- data for evaluation of derivatives ----
get_inner_deriv_data <- function(xdat, model_data, betasim, term_names, deriv_terms, model_terms) {

  colnames(xdat) <- term_names$covariates

  betasim <- data.frame(betasim)
  names(betasim) <- term_names$betas
  n_betasim <- nrow(betasim)

  # add data of constitutive covariates if necessary
  # put this into specific function?
  term_labels <- attr(model_terms, "term.labels")
  int <- term_labels[grepl(":", term_labels)]
  int_terms <- unlist(strsplit(int, ":", fixed = TRUE))

  check_terms <- int_terms %in% term_labels == FALSE

  if (any(check_terms)) {
    addvars <- int_terms[check_terms]
    xdat <- cbind(xdat, as.matrix(model_data[, addvars, drop = FALSE]))
  }

  # keep only data required to compute derivatives (reduce memory use)
  sel_xdat <- colnames(xdat)[colnames(xdat) %in% deriv_terms]
  sel_betasim <- names(betasim)[names(betasim) %in% deriv_terms]

  xdat <- xdat[, sel_xdat, drop = FALSE]
  betasim <- betasim[, sel_betasim, drop = FALSE]

  # set up list of dataframes needed to compute derivatives
  inner_deriv_data <- vector(mode = "list", length = n_betasim)
  for (i in 1:n_betasim) {
    beta <- betasim[i, , drop = FALSE]
    inner_deriv_data[[i]] <- cbind(xdat, beta)
  }
  return(inner_deriv_data)
}

## ---- function to evaluate inner derivatives ----
eval_inner_derivs <- function(d, n_sim, inner_deriv_data, dydx_vars) {
  n_d <- length(d)
  n_data <- nrow(inner_deriv_data[[1]])

  dydx <- vector(mode = "list", length = n_d)
  # loop over dydx_vars
  for (i in 1:n_d) {
    grad <- matrix(NA, nrow = n_data, ncol = n_sim)
    # loop over simulations. Evaluate partial derivatives
    for (j in 1:n_sim) {
      grad[,j] <- eval(d[[i]], inner_deriv_data[[j]])
    }
    dydx[[i]] <- grad
  }
  names(dydx) <- dydx_vars
  return(dydx)
}


