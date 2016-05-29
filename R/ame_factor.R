ame_factor <- function(fac_vars, coefnames, betasim, xdat, model_data, model_family, model_link) {

  # what does the function do?
  # 1. identify variables relevant for computation of first differences:
  #    - factor variables
  #    - continuous variables in factor*continuous interaction terms (how about factor * factor interactions?)
  # 2. set scenarios (base line, first factor level, second factor level, etc.)
  #    - set (constitutive) dummy factor terms in model matrix to zero/one as required
  #    - extract constitutive continuous variables of factor*continuous interaction terms from model frame
  #    - recompute factor dummy*continuous interaction terms in model matrix
  # 3. compute first differences:
  #    - compute linear predictors for each scenario
  #    - apply link function to linear predictor
  #    - get difference between base factor level and other factor levels

  # 1. identify relevant variables/model terms
  n_fac_vars <- length(fac_vars)

  # get all terms containing a factor variable
  contains_factor <- lapply(fac_vars, function(fv) coefnames[grep(fv, coefnames)])

  # get factor levels
  factor_levels <- vector(mode = "list", length = n_fac_vars)
  for (i in 1:n_fac_vars) {
    fl <- contrasts(model_data[, fac_vars[i]])
    factor_levels[[i]] <- paste0(fac_vars[i], rownames(fl))
  }

  # get constitutive factor terms
  factor_terms <- vector(mode = "list", length = n_fac_vars)
  for (i in 1:n_fac_vars) {
    factor_terms[[i]] <- contains_factor[[i]][grep(":", contains_factor[[i]], fixed = TRUE, invert = TRUE)]
  }

  # get interaction terms
  int_terms <- vector(mode = "list", length = n_fac_vars)
  for (i in 1:n_fac_vars) {
    int_terms[[i]] <- contains_factor[[i]][grep(":", contains_factor[[i]], fixed = TRUE)]
  }

  # get continuous variables of factor * continuous interaction terms
  cont_terms <- vector(mode = "list", length = n_fac_vars)
  for (i in 1:n_fac_vars) {
    c_t <- unique(unlist(strsplit(contains_factor[[i]], ":", fixed = TRUE)))
    cont_terms[[i]] <- c_t[c_t %in% factor_levels[[i]] == FALSE]
  }

  # get base levels
  base_levels <- vector(mode = "list", length = n_fac_vars)
  for (i in 1:n_fac_vars) {
    base_levels[[i]] <- factor_levels[[i]][factor_levels[[i]] %in% factor_terms[[i]] == FALSE]
  }

  # 2. set scenarios
  factor_xdat <- vector(mode = "list", length = n_fac_vars)
  for (i in 1:n_fac_vars) {

    f_l <- factor_levels[[i]]
    c_f <- contains_factor[[i]]
    f_t <- factor_terms[[i]]
    c_t <- cont_terms[[i]]
    i_t <- int_terms[[i]]
    b_l <- base_levels[[i]]

    n_fac_levels <- length(f_l)
    all_levels <- c(b_l, f_t)

    # set scenarios for each factor variable
    xdat_list <- vector(mode = "list", length = n_fac_levels)
    for (j in 1:n_fac_levels) {

      # get model matrix
      xdat_list[[j]] <- xdat

      # complementary factor levels
      c_fl <- c_f[grep(all_levels[j], c_f, fixed = TRUE, invert = TRUE)]

      # set complementary factor levels to zero
      xdat_list[[j]][, c_fl] <- 0

      # set main factor level to one; recompute interaction terms
      m_fl <- all_levels[j]
      if (j > 1) {
        xdat_list[[j]][, m_fl] <- 1

        if (length(i_t) > 0) {
          # get interaction terms of interest
          fxc0 <- i_t[grep(m_fl, i_t)]
          fxc1 <- gsub(paste0(m_fl, ":"), "", fxc0)
          fxc2 <- gsub(":", "*", fxc1)

          # recompute interaction terms
          fxc_call <- lapply(fxc2, function(x) parse(text = x))
          x <- sapply(fxc_call, eval, model_data)
          xdat_list[[j]][, fxc0] <- x
        }
      }
    }
    factor_xdat[[i]] <- xdat_list
  }

  # 3. compute first differences
  # linear predictor mu
  mu <- vector(mode = "list", length = n_fac_vars)
  for (i in 1:n_fac_vars) {
    mu[[i]] <- lapply(factor_xdat[[i]], function(x) x %*% t(betasim))
  }

  # expected values
  if (model_family == "gaussian" & model_link == "identity") {

    f_mu <- mu

  } else if (model_family == "binomial" & model_link == "logit") {

    f_mu <- vector(mode = "list", length = n_fac_vars)
    for (i in 1:n_fac_vars) {
      f_mu[[i]] <- lapply(mu[[i]], plogis)
    }

  } else if (model_family == "binomial" & model_link == "probit") {

    f_mu <- vector(mode = "list", length = n_fac_vars)
    for (i in 1:n_fac_vars) {
      f_mu[[i]] <- lapply(mu[[i]], pnorm)
    }

  } else {
    stop("link function not specified")
  }

  # compute first differences
  fd_list <- vector(mode = "list", length = n_fac_vars)
  for (i in 1:n_fac_vars) {

    # get f_mu for each factor variable
    fmu <- f_mu[[i]]
    l_fmu <- length(fmu)

    # compute first differences
    first_difference <- vector(mode = "list", length = l_fmu - 1)
    names(first_difference) <- factor_terms[[i]]
    for (j in 2:l_fmu) {
      first_difference[[j-1]] <- colMeans(fmu[[j]] - fmu[[1]])
    }

    fd_list[[i]] <- first_difference
  }

  # simplify fd_list
  fd_list <- unlist(fd_list, recursive = FALSE)

  # set attributes
  attr(fd_list, "variable_type") <- "factor"
  attr(fd_list, "variable_names") <- fac_vars

  return(fd_list)

}
