## ---- documentation of ame ----
#' Average Marginal Effects
#'
#' Compute average marginal effects for generalized linear and generalized additive models using Krinsky-Robb Monte-Carlo simulation.
#'
#' @author Julia Partheymueller, Alexander Staudt (\email{astaudt@@mail.uni-mannheim.de})
#'
#' @param model an object of class \code{glm}.
#'
#' @param cont_vars a character vector containing the names of continuous model variables.
#'
#' @param fac_vars a character vector with factor model variables.
#'
#' @param nsim number of simulations. Default is \code{100}.
#'
#' @param seed Seed value for the random number generator. Default value is \code{1}.
#'
#' @return \code{ame} returns an object of class \code{"ame"}. The functions \code{summary} (\code{summary.ame}) and \code{print} (\code{print.summary.ame}) are used to display a summary of the computed average marginal effects.
#'
#' An object of class \code{"ame"} is a list that contains at least the following components:
#'
#' \item{\code{continuous}}{a named list containing average marginal effects for each continuous variable. If no variable is specified, defaults to \code{NULL}. The list has attributes \code{variable_type} and \code{variable_names}.}
#' \item{factor}{a named list containing first differences for each factor variable. If no factor is specified, defaults to \code{NULL}. The list has attributes \code{variable_type} and \code{variable_names}.}
#'
#' Furthermore, objects of class \code{"ame"} store the model family and link in the attributes \code{model_family}, \code{model_link}, and the random number generator seed and number of draws from the multivariate normal in the attributes \code{seed} and \code{n_sim}.
#'
#' @references
#'
#' Krinsky I, Robb AL (1991). "Three Methods for Calculating the Statistical Properties of Elasticities: A Comparison." \emph{Empirical Economics} 16(2), 199--209.
#'
#' Brambor T, Clark WR, Golder M (2006). "Understanding Interaction Models: Improving Empirical Analyses." \emph{Political Analysis} 14(1), 63--82. url: \url{http://localgov.fsu.edu/readings_papers/Research Methods/Brambor_et_al_Multipolicative_Interactions.pdf}
#'
#' @seealso
#' King G, Tomz M, Wittenberg J (2000). "Making the Most of Statistical Analyses: Improving Interpretation and Presentation." \emph{American Journal of Political Science} 44(2), 347--361. url: \url{http://www.polmeth.wustl.edu/files/polmeth/king98f.pdf}.
#'
#' @examples
#' # generate some data (factor, dummy, continuous)
#' set.seed(100)
#'
#' y <- rbinom(100, size = 1, prob = 0.5)
#' mcat <- as.factor(rpois(100, 3))
#' dummy <- as.factor(rbinom(100, size = 1, prob = 0.5))
#' cont <- runif(100, -10, 10)
#' data <- data.frame(y, mcat, dummy, cont)
#'
#' # run glm model
#' glm_fit <- glm(y ~ mcat + dummy * cont, data = data, family = binomial(link = "logit"))
#'
#' summary(glm_fit)
#'
#' # compute average marginal effects
#' glm_ame <- ame(glm_fit, cont_vars = c("cont"), fac_vars = c("mcat", "dummy"), nsim = 1000)
#'
#' summary(glm_ame)
#'
#' # linear model using glm
#' lm1 <- glm(mpg ~ cyl + hp + wt, data = mtcars, family = gaussian)
#' lm2 <- glm(mpg ~ cyl * hp + wt, data = mtcars, family = gaussian)
#'
#' summary(lm1)
#' summary(lm2)
#'
#' # compute average marginal effects
#' lm1_ame <- ame(lm1, cont_vars = c("cyl", "hp", "wt"), nsim = 1000, seed = 99)
#' lm2_ame <- ame(lm2, cont_vars = c("cyl", "hp", "wt"), nsim = 1000, seed = 99)
#'
#' # display results
#' summary(lm1_ame)
#' summary(lm2_ame)
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats D coef contrasts dlogis dnorm model.frame model.matrix plogis pnorm printCoefmat reformulate sd terms vcov
#' @export

## ---- ame function ----
ame <- function(model, cont_vars, fac_vars, nsim = 100, seed = 1) {
  if (missing(model)) {
    stop("Please specify your model.")
  }
  if (missing(cont_vars) & missing(fac_vars)) {
    stop("Please specify your variables of interest.")
  }

  # check if object has class "gam"
  if ("gam" %in% class(model) | "glm" %in% class(model)) {

    # get distribution and link function
    model_family <- model$family$family
    model_link <- model$family$link

    # Step 1: get coefficients and variance-covariance of the model
    coefs <- coef(model)
    vc <- vcov(model)

    # coefficient names
    coefnames <- names(coefs)

    # Step 2: simulate coefficients from multivariate normal
    set.seed(seed)
    betasim <- mvrnorm(nsim, mu = coefs, Sigma = vc)

    # Step 3: get model matrix and model terms
    # does not deal with I() statements yet
    xdat <- model.matrix(model)
    model_data <- model.frame(model)

    model_terms <- terms(model)


    # Calculate marginal effects and first differences
    if (!missing(cont_vars)) {
      results_continuous <- ame_continuous(cont_vars = cont_vars,
                                           n_sim = nsim,
                                           coefnames = coefnames,
                                           betasim = betasim,
                                           xdat = xdat,
                                           model_data = model_data,
                                           model_family = model_family,
                                           model_link = model_link,
                                           model_terms = model_terms)
    }
    else {
      results_continuous <- NULL
    }
    if (!missing(fac_vars)) {
      results_factor <- ame_factor(fac_vars = fac_vars,
                                   coefnames = coefnames,
                                   betasim = betasim,
                                   xdat = xdat,
                                   model_data = model_data,
                                   model_family = model_family,
                                   model_link = model_link)
    }
    else {
      results_factor <- NULL
    }

    results <- list(continuous = results_continuous, factor = results_factor)
    attr(results, "model_family") <- model_family
    attr(results, "model_link") <- model_link
    attr(results, "n_sim") <- nsim
    attr(results, "seed")  <- seed

    class(results) <- "ame"

    return(results)

  } else {
    stop("Model class not supported.")
  }

}

