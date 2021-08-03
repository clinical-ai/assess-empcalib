## ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##
## Functions to perform modelling
##
## ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###




#' Generates additive formula for non treatment covariates.
#' E.g., \code{C1 + C2 + C3}
#'
#' @param data Data to extract covariates from
#' @param covar_prefix Prefix of covariate columns
#'
#' @return Character vector containing the model formula. E.g., \code{C1 + C2 + C3}
#'
#' @export
gen_additive_non_trt_covars_formula <- function(data = NULL, covar_prefix = "C") {
    # Additive terms for the covariates
    baseline_covars <- data %>% select(starts_with(covar_prefix)) %>% colnames()
    rhs_formula <- paste0(baseline_covars, collapse = "+")
    return(rhs_formula)
}

#' Generates propensity scoring formula with quadratic term.
#' E.g., \code{Z ~ X1 + X1^2 + X2 + X2^2}
#'
#' @param trt Treatment variable
#' @param covar_prefix Prefix of baseline (non-treatment) covariates
#' @param total_num_covars Total number of covariates
#'
#' @return Character vector of quadratic propensity score formula.
gen_pscore_model_formula_proper_quadratic <- function(trt = "Z", covar_prefix = "C", total_num_covars = 1L) {
}

#' Propensity score model formula with missing confounders
#'
#' @export
gen_pscore_model_formula_ignoreconfounder <- function(trt = "Z", covar_prefix = "C",
                                                      total_num_covars = 2L, covar_incre_by = 2L) {
    rhs_formula <- gen_formula_rhs_with_gaps(covar_prefix, total_num_covars, covar_incre_by)
    return(paste(trt, rhs_formula, sep = " ~ "))
}


