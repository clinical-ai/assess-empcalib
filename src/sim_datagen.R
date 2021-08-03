library(dplyr)
library(tibble)
library(glue)
library(here)

#' Creates a data table that contains data generation specifications for the confounders.
#' Part of v2 architecture.
#'
#' @param n Denotes the number of confounders to create.
#' @param datagen_fn Name of function to generate values.
#' @param genfn_args Arugments of data generation function specified in \code{datagen_fn}
#' @param confounder_name_prefix
create_confounder_spec <- function(n = 0L,
                                   datagen_fn = "rnorm",
                                   genfn_args = list(mean = 0),
                                   confounder_name_prefix = "X") {

    if (!require(purrr)) { library(purrr) }

    num_cols <- n

    # Constructs specification table to construct confounders
    confounder_spec <- seq.int(from = 1L, to = num_cols) %>% purrr::map_dfr(.f = ~{
        return(list(
            confounder_num = .x,
            genfn = datagen_fn,
            genfn_args = list(genfn_args),
            confounder_name = paste0(confounder_name_prefix, .x))
        )
    })

    return(confounder_spec)

} # End gen_confounder_values()

add_one_confounder_spec <- function(spec = NULL,
    datagen_fn = "rnorm", genfn_args = list(mean = 0),
    covar_name = "X1") {

    if (!require(tibble)) { library(tibble) }

    if (is.null(spec)) { return(NULL) }

    new_id <- (spec %>% pull(confounder_num) %>% max()) + 1L
    new_spec_row <- create_confounder_spec(n = 1L, datagen_fn = datagen_fn,
        genfn_args = genfn_args, confounder_name_prefix = covar_name)

    new_spec_row <- new_spec_row %>% mutate(confounder_num = c(new_id))
    new_spec_row <- new_spec_row %>% mutate(confounder_name = c(covar_name))

    spec_to_rt <- spec %>% tibble::add_row(new_spec_row)
} # End `add_one_confounder_spec()`

gen_confounders_from_spec <- function(n_obs = 0L, confounder_spec = NULL) {

    if (!require(purrr)) { library(purrr) }

    confounders_mat <- confounder_spec %>% purrr::pmap_dfc(.f = ~{
        colnum <- ..1
        datagen_fn <- ..2
        datagen_fn_args <- ..3
        confounder_name <- ..4
        n_obs <- ..5
        datagen_fn_args[['n']] <- n_obs

        vals <- do.call(what = datagen_fn, args = datagen_fn_args)

        rt_list <- list()
        rt_list[[confounder_name]] = vals
        return(rt_list)

    }, n = n_obs)

    return(confounders_mat)

} # End gen_confounders_from_spec()

#' Generate values of regression coefficients from confounder specification table.
#' Part of v2 arch.
gen_trt_model_regression_coefs <- function(confounder_spec = NULL,
    coef_range = c(min_or = 1.1, max_or = 4.0),
    coef_to_sample = NULL,
    # coef_to_sample = c(log(1.5), log(2.0), log(4.0), log(6.0)),
    incl_intercept = TRUE) {

    if (!require(vctrs)) { library(vctrs) }
    if (is.null(confounder_spec)) { return(c(0.0)) }

    n_coefs <- ifelse(TRUE == incl_intercept,
        yes = nrow(confounder_spec) + 1L, no = nrow(confounder_spec))

    if (is.null(coef_to_sample) & !is.null(coef_range)) {
        trt_model_coef_vec <- log( runif(n_coefs,
                                         min = coef_range[['min_or']], max = coef_range[['max_or']]))

    } else if (!is.null(coef_to_sample) & is.null(coef_range)) {
        # From ?sample() this is to take care of case where there is only
        # one coefficient to sample.
        resample <- function(x, ...) x[sample.int(length(x), ...)]
        trt_model_coef_vec <- resample(coef_to_sample,
            size = n_coefs,
            replace = TRUE)
    }

    coef_names <- c()
    if (TRUE == incl_intercept) {
        coef_names <- c("INTR", confounder_spec %>% pull(confounder_name))
    } else {
        coef_names <- confounder_spec %>% pull(confounder_name)
    }

    names(trt_model_coef_vec) <- coef_names
    return(trt_model_coef_vec)
}


#' Generates a matrix containing regression coefficients of confounder to
#' outcomes.
#' Part of v2 architecture
#'
#' @param total_num_negctrls Number of negative control outcomes
#' @param num_not_pri_outcome_confounders Integer that indicates the number of
#'   confounders not associated with primary outcome (outcome of interest).
#' @param trt_eff_logor_outcomeinterest denotes the regression coefficient
#'   of treatment to outcome of interest, this is in \emph{log odds ratio}.
#' @param trt_eff_range_or Vector with the range of regression coefficient
#'   of treatment to outcome of interest in \emph{odds ratio}, specified by
#'   two named elements \code{min_or} and \code{max_or}.
#' @param coef_prefix Prefix character for confounders
#' @param coef_prefix Variable name for treatment value
gen_outcome_model_coef_template <- function(
    total_num_negctrls = 0,
    num_not_pri_outcome_confounders = 0L,
    trt_eff_logor_outcomeinterest = NA_real_,
    trt_eff_range_or = NULL,
    coef_prefix = "b", trt_valname = "z",
    outcome_id_colname = "outcome_id") {

    if (!require(tibble)) { library(tibble) }
    if (!require(purrr))  { library(purrr) }

    num_outcomes <- total_num_negctrls + 1L # Plus 1 for outcome of interest

    # Negative controls have no treatment effect
    if (!is.na(trt_eff_logor_outcomeinterest)) {
        trt_coefs <- c(trt_eff_logor_outcomeinterest, rep(0.0, times = total_num_negctrls))
    } else if (!is.null(trt_eff_range_or)) {
        trt_coefs <- c( log(runif(1 , min = trt_eff_range_or[['min_or']],
                                      max = trt_eff_range_or[['max_or']])),
                        rep(0.0, times = total_num_negctrls))
    } else {
        trt_coefs <- NA_real_
    }

    # Initial table of coefficients; first row concerns primary outcome.
    outcome_coef_tbl <- tibble(
        !!outcome_id_colname := seq.int(1L, to = num_outcomes),
        !!paste0(coef_prefix, 0) := rep(NA_real_, times = num_outcomes),
        !!paste0(coef_prefix, trt_valname) := trt_coefs,
        !!paste0(coef_prefix, 1) := rep(NA_real_, times = num_outcomes)
    )

    notpriout_conf_tbl <- seq.int(2L, to = num_not_pri_outcome_confounders + 1L) %>%
        purrr::map_dfc( .f = ~ {
            coef_name <- paste0(coef_prefix, .x)
            list(coef_name = rep(NA_real_, times = num_outcomes))
        })
    names(notpriout_conf_tbl) <- paste0(coef_prefix,
        seq.int(2L, to = num_not_pri_outcome_confounders + 1L))

    outcome_coef_tbl <- cbind(outcome_coef_tbl, notpriout_conf_tbl)
    return(outcome_coef_tbl)

} # End gen_outcome_model_coef_template()

add_one_outcome_model_coef <- function(outcome_model_coefs = NULL,
    coef_name = "U", coefs_vals = c()) {

    if (!require(dplyr)) { library(dplyr) }

    if (is.null(outcome_model_coefs)) { return(NULL) }

    coef_vals_to_add <- rep(NA_real_, times = nrow(outcome_model_coefs))
    if ( (length(coefs_vals) != 0) & (FALSE == is.null(coefs_vals))) {
        if (length(coefs_vals) != nrow(outcome_model_coefs)) {
            stop("!! ERROR: add_one_outcome_model_coef() Incorrect length of supplied coef values vector.")
        } else {
            coef_vals_to_add <- coefs_vals
        }
    }

    outcome_model_coefs <- outcome_model_coefs %>% mutate(
        !!coef_name := coef_vals_to_add)

    return(outcome_model_coefs)

} # End add_one_outcome_model_coef()

#' Generate values for coefficient vector. Arch v2.
#' Assumes the first coefficient is for the outcome of interest. Therefore it's value is
#' always set by default (use \code{outcome_of_interest_always_set} to change this
#' behaviour).
#' Proportion to set is then applicable to the coefficients for the negative controls by default.
set_coef_vals <- function(initial_coef_vec = c(), prop_negctrl_to_set = 1.0,
    val_range = c(),
    val = NA_real_,
    outcome_of_interest_always_set = TRUE,
    coef_noeff_val = NA_real_) {

    if (!require(dplyr)) { library(dplyr) }

    if (0 == length(initial_coef_vec)) { return(initial_coef_vec) }

    resample <- function(x, ...) x[sample.int(length(x), ...)]

    num_negctrl_indicies_to_set <- ceiling( (length(initial_coef_vec) - 1) * prop_negctrl_to_set) %>% as.integer()
    rand_negctrl_indicies_to_set <- resample(seq.int(from = 2L, to = length(initial_coef_vec)),
            size = num_negctrl_indicies_to_set)

    coef_to_rt <- initial_coef_vec
    if (!is.na(val)) {
        coef_to_rt[rand_negctrl_indicies_to_set] <- val
    } else {
        if (0 == length(val_range)) {
            stop("!! ERROR: set_coef_vals() No coefficient range to sample from")
        }
        # vals_to_set <- resample(val_range, size = num_negctrl_indicies_to_set, replace = TRUE)
        vals_to_set <- log( runif(n = num_negctrl_indicies_to_set,
                                  min = val_range[['min_or']], max = val_range[['max_or']] ))
        coef_to_rt[rand_negctrl_indicies_to_set] <- vals_to_set
    }

    outcome_coef_not_set <- FALSE
    if (is.na(coef_noeff_val)) {
        outcome_coef_not_set <- is.na(coef_to_rt[[1]])
    } else {
        outcome_coef_not_set <- (coef_to_rt[[1]] == coef_noeff_val)
    }

    if (outcome_of_interest_always_set & outcome_coef_not_set) {
        # Coefficient for outcome of interest needs to be set.
        if (!is.na(val)) {
            coef_to_rt[[1]] <- val
        } else {
             # coef_to_rt[[1]] <- resample(val_range, size = 1)
            coef_to_rt[[1]] <- log( runif(n = 1, min = val_range[['min_or']], max = val_range[['max_or']] ))
        }
    }

    # Sets the remaining indices to zero.
    indicies_with_vals <- c(1L, rand_negctrl_indicies_to_set)
    if (sum (is.na(coef_to_rt[-indicies_with_vals]) ) > 0) {
        coef_to_rt[-indicies_with_vals] <- 0.0
    }

    return(coef_to_rt)

} # End gen_potential_confounder_outcome_coef_vec()

#' Create negative control data by combining the vector of outcomes and covariates.
#' Part of equi-confounding, aka v2 architecture.
create_negctrl_data <- function(outcome_id = NA_integer_, negctrl_outcomes_mat = NULL,
                                covars_data = NULL,
                                outcome_var_name = "Y") {

    if (is.na(outcome_id) | is.null(negctrl_outcomes_mat) | is.null(covars_data)) {
        error("!! ERROR: create_negctrl_data() One of required data is NULL or NA.")
    }

    if (!require(dplyr)) { library(dplyr) }

    # Assumes `negctrl_outcomes_mat` is a matrix with dimension 'n' by 'num_negctrl + 1',
    # 'n' is the number of observations, and the first column is the outcome of interest.
    negctrl_outcome_vec <- negctrl_outcomes_mat[, outcome_id]

    # Assumes 'trt_data' is a data.frame like structure that contains all the covariates
    # appropriate for negative control—i.e., columns such as potential unmeasured confounders
    # are removed.
    negctrl_data <- covars_data %>% mutate(!!outcome_var_name := negctrl_outcome_vec) %>%
        relocate(.data[[outcome_var_name]])

    return(list(outcome_id = outcome_id,
                negctrl_id = outcome_id - 1L,
                negctrl_data = list(negctrl_data)))

} # End create_negctrl_data()


#' Obtains non-treatment covariate column names from data.
#'
#' @param coef_tbl
#' @param coef_prefix
#' @return Vector of characters
get_nontrt_covarnames <- function(coef_tbl = NULL, coef_prefix = "b") {
    covar_names <- grep(glue("^{coef_prefix}[0-9]+$"),
                        x = names(coef_tbl), value = TRUE)

} # End of get_nontrt_covarnames()

#' Sets all \emph{non-treatment} co-efficients across all outcomes.
#' E.g., when in use with \code{use_outcome_model_coef_template_all_sims} flag.
#' Part of v2 architecture.
#'
#' @param coef_template_tbl Coefficient template table. Assumes each row
#'   is an outcome, each column is a covariate. Covariates are named
#'   as \code{coef_prefix<number>} format, where \code{number} is an increasing
#'   integer. E.g., \code{b1, b2, ...}
#' @param coef_prefix
#' @param coef_val_range Range of values to sample new coefficient from.
#'   Requires \code{min_or} and \code{max_or}.
set_nontrt_varcoefs <- function(coef_template_tbl = NULL,
                                coef_prefix = "b",
                                coef_val_range = c(min_or = 1.1, max_or = 4.0)) {

    covar_names <- get_nontrt_covarnames(coef_template_tbl, coef_prefix)

    for (measured_conf_colname in covar_names) {

        current_coef_vec <- coef_template_tbl %>% pull(.data[[measured_conf_colname]])
        beta_x_vec <- set_coef_vals(initial_coef_vec = current_coef_vec,
                                        prop_negctrl_to_set = 1.0,
                                        val_range = coef_val_range)

        coef_template_tbl <- coef_template_tbl %>%
            mutate(!!measured_conf_colname := beta_x_vec)
    }
    return(coef_template_tbl)

} # End set_nontrt_varcoefs()

#' Sets all \emph{non-treatment} co-efficients across all outcomes
#' \emph{to the same value}.
#' E.g., when in use with \code{use_outcome_model_coef_template_all_sims} flag.
#' Part of v2 architecture.
#'
#' @param coef_template_tbl Coefficient template table. Assumes each row
#'   is an outcome, each column is a covariate. Covariates are named
#'   as \code{coef_prefix<number>} format, where \code{number} is an increasing
#'   integer. E.g., \code{b1, b2, ...}
#' @param coef_prefix
#' @param coef_val Value to set.
#'   Requires \code{min_or} and \code{max_or}.
set_nontrt_varcoefs_sameval <- function(coef_template_tbl = NULL,
                                        coef_prefix = "b",
                                        coef_val = log(1.0)) {

    covar_names <- get_nontrt_covarnames(coef_template_tbl, coef_prefix)

    for (measured_conf_colname in covar_names) {
        current_coef_vec <- coef_template_tbl %>% pull(.data[[measured_conf_colname]])
        beta_x_vec <- set_coef_vals(initial_coef_vec = current_coef_vec,
                                    prop_negctrl_to_set = 1.0,
                                    val = coef_val)

        coef_template_tbl <- coef_template_tbl %>%
            mutate(!!measured_conf_colname := beta_x_vec)
    }
    return(coef_template_tbl)

} # End set_nontrt_varcoefs_sameval()

#' Generates positive control data, using regression estimates of
#' negative control. Part of new architecture to handle equi-confounding. v2 arch
gen_posctrl_data_v2 <- function(target_trt_eff = NA_real_, data = NULL, coef_template_vec = NULL,
                             coef_prefix = "c",
                             trt_varname = "Z", outcome_varname = "Y") {

    if (!require(dplyr)) { library(dplyr) }
    if (!require(glue))  { library(glue) }

    if (target_trt_eff < 0) {
        warning("!! WARNING: gen_posctrl_data() Target treatment effect is less than zero.")
    }
    if (is.null(coef_template_vec)) {
        error("!! ERROR: gen_posctrl_data() Coefficient template not supplied.")
    }
    if (is.null(data)) {
        error("!! ERROR: gen_posctrl_data() Confounder data not supplied.")
    }

    print(glue("> DEBUG: gen_posctrl_data() Generating positive control at level {signif(target_trt_eff, 3)}"))

    # Assumes the order of columns in `data` and `coef_template` are the same.
    # E.g.,
    # Z, X1, X2, ... , U1 , etc (i.e., treatment, measured confounders, extra terms)

    # Creates design matrix
    design_mat <- data %>% select(-.data[[outcome_varname]]) %>% # Removes neg ctrl outcome column
        mutate(INTR = rep(1.0, times = nrow(data))) %>%          # Adds intercept term
        relocate(INTR) %>%
        as.matrix()

    # Sets the correct level for treatment
    coef_template_vec[[paste0(coef_prefix, trt_varname)]] <- target_trt_eff

    # Outcome vector
    outcome_model_eta <- design_mat %*% coef_template_vec
    pr_y1 <- 1 / (1 + exp(-1 * outcome_model_eta))
    y_vec <- rbinom(nrow(data), size = 1L, prob = pr_y1)

    postctrl_data <- data %>% mutate(!!outcome_varname := y_vec)

    return(list(trt_eff = target_trt_eff, postctrl_data = list(postctrl_data)))

} # End gen_posctrl_data()

#' Returns the range of odds-ratios for coefficients
#'
#' @return Vector with minimum and maximum \emph{odds ratio} stored in
#'     \code{min_or} and \code{max_or} respectively.
get_coef_range_or <- function() {
    coef_range <- c(min_or = 0.5, max_or = 2)
    return(coef_range)
}
