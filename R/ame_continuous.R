ame_continuous <- function(cont_vars, n_sim, coefnames, betasim, xdat, model_data, model_family, model_link, model_terms){

  # 1. get linear predictor formula (mu)
  # coefficient names
  term_names <- get_term_names(coefnames)

  # get linear function term (mu)
  mu <- paste(term_names$betas, term_names$covariates, sep = " * ", collapse = " + ")

  # 2. obtain partial inner derivatives of mu
  d <- vector(mode = "list", length = length(cont_vars))
  for (i in 1:length(cont_vars)){
    d[[i]] <- D(reformulate(mu)[[2]], cont_vars[i])
  }

  deriv_terms <- get_deriv_terms(d, cont_vars)

  # 3. compute partial inner derivatives (inner dydx)
  # get data for function evaluation of analytic derivatives
  xdat_beta <- get_inner_deriv_data(xdat, model_data, betasim, term_names, deriv_terms, model_terms)

  # evaluate analytical inner derivatives
  inner_dydx <- eval_inner_derivs(d, n_sim, xdat_beta, cont_vars)

  # 4. compute outer derivative and marginal effects for specified variables
  if (model_family == "gaussian" & model_link == "identity") {

    l <- 1
    # marginal effect averaged over all observations for each draw of betasim
    me_coefsmean <- lapply(inner_dydx, function(x) colMeans(x * l))

  } else if (model_family == "binomial" & model_link == "logit") {

    l <- apply(betasim, 1, function(x) dlogis(xdat %*% x))
    me_coefsmean <- lapply(inner_dydx, function(x) colMeans(x * l))

  } else if (model_family == "binomial" & model_link == "probit") {

    l <- apply(betasim, 1, function(x) dnorm(xdat %*% x))
    me_coefsmean <- lapply(inner_dydx, function(x) colMeans(x * l))

  } else {

    stop("link function not specified")

  }

  # set attributes
  attr(me_coefsmean, "variable_type") <- "continuous"
  attr(me_coefsmean, "variable_names") <- cont_vars

  return(me_coefsmean)

}
