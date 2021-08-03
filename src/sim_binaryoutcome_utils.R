## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##
## Supporting functions for the binary outcome model
##
## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#' @export
get_outcome_of_interest_study_num <- function() { return(1L) }

#' @export
get_studyid_colname <- function() { return("outcomeid") }

#' @export
get_intercept_term_colname <- function() { return("INTR") }

#' @export
get_extra_confounder_colname <- function() { return("U") }

#' @export
get_confounder_colname_prefix <- function() { return("X") }

#' @export
get_outcomedatatype_colname_prefix <- function() { return("type") }

#' @export
get_subjid_colname_prefix <- function() { return("ID") }

#' @export
get_outofint_confounder_colname <- function() { return( paste0(get_confounder_colname_prefix(), "1")) }

#' @export
get_trt_colname <- function() { return("Z") }

#' @export
get_outcome_colname <- function() { return("Y") }

#' @export
get_quadterm_colname <- function() { return("X1Sqr") }

#' @export
get_x1x2intr_colname <- function() { return("X1X2") }

#' @export
get_x1_zintr_colname <- function() { return("X1Z") }

#' @export
get_trtpos_deltavar_colname <- function() { return("XTRT") }

#' @export
get_merror_confounder_colname <- function() { return("XERR") }

#' Returns the strings associated with different types of calibrations
get_syserr_model_calibration_type <- function() {
    return(list(nomod = "nomod",
                negctrlonlysyserr = "negctrlonlysyserr",
                nullsyserr = "nullsyserr",
                adjtruthposctrlsyserr = "adjtruthposctrlsyserr"))
}

gen_negctrl <- function(x1_vec = c(), x2_vec = c(), x3_vec = c(), u_vec = c(), z_vec = c(),
                        beta_vec = c(),
                        incl_covar_u = FALSE,
                        incl_x1sqr = FALSE,
                        incl_x1x2intr = FALSE,
                        incl_x1zintr = FALSE,
                        outcome_generator = gen_outcome_vec) {

    b0 <- beta_vec[[1]] # b0_s
    b1 <- beta_vec[[2]] # b1_s
    b2 <- beta_vec[[3]] # b2_s
    b3 <- beta_vec[[4]] # b3_s

    bz <- log(1.0) # odds ratio of 1, p = 0.5

    b1sqr <- ifelse(gen_quad_term_in_negctrl(), yes = beta_vec[[5]], no = NA_real_)

    bu <- ifelse(gen_extra_confounder_in_negctrl(), yes = beta_vec[[6]], no = NA_real_)

    bx1x2 <- ifelse(gen_x1x2_interaction_in_negctrl(), yes = beta_vec[[7]], no = NA_real_)
    bx1z  <- ifelse(gen_x1z_interaction_in_study(), yes = beta_vec[[8]], no = NA_real_)

    Y_data <- gen_study_data(is_negctrl = TRUE,
                             x1_vec = x1_vec, x2_vec = x2_vec, x3_vec = x3_vec, u_vec = u_vec, z_vec = z_vec,
                             intercept_coef = b0,
                             b1 = b1, b2 = b2, b3 = b3, b_u = bu, b_z = bz,
                             incl_covar_u = incl_covar_u,
                             incl_x1sqr = incl_x1sqr, b_x1sqr = b1sqr,
                             incl_x1x2intr = incl_x1x2intr, b_x1x2 = bx1x2,
                             incl_x1zintr = incl_x1zintr, b_x1z = bx1z,
                             outcome_generator = outcome_generator)
    return(Y_data)

} # End `gen_negctrl()`

#' Generates data for positive control at a specific target odds ratio
#'
#' @param negctrl_data Data of current negative control. Can contain columns suhc as
#'   subject ID, and negative control outcome.
#' @param est_coef_negctrl_data_vec Named vector of estimated regression coefficients from
#'   negative control data specified in \code{negctrl_data} argument.
#'
#' @return \link{\code{tibble}} of positive control data, including subject IDs,
#'   treatment and new outcome.
#'
#' @export
gen_posctrl_data <- function(negctrl_data = NULL,
                             est_coef_negctrl_data_vec = c(NA_real_),
                              target_odds_ratio = 2.0,
                              outcome_generator = NULL) {

    est_coef_negctrl_data_vec[[get_trt_colname()]] <- log(target_odds_ratio)

    # Removes columns that are not part of design matrix.
    #print("DEBUG gen_posctrl_data(): Removing columns")

    cols_to_remove_negctrl <- vec_c(get_outcomedatatype_colname_prefix(),
                                    get_subjid_colname_prefix(),
                                    get_outcome_colname())
    current_negctrl_design_matrix <- negctrl_data %>%
        select(-any_of(cols_to_remove_negctrl)) %>% as.matrix()

    # Constructs design matrix.
    #print("DEBUG gen_posctrl_data(): Constructing design matrix")
    intercept_vec <- rep(1L, times = nrow(current_negctrl_design_matrix))
    current_negctrl_design_matrix <- cbind(INTR = intercept_vec,
                                           current_negctrl_design_matrix)

    # Obtains linear combination of covariates.
    liner_combo_vec <- current_negctrl_design_matrix %*% est_coef_negctrl_data_vec
    # Generates outcomes using linear combination of covariates.
    #print("DEBUG gen_posctrl_data(): Generating positive control outcomes")
    postctrl_outcomes_vec <- do.call(what = outcome_generator,
                                     args = list(n_obs = nrow(current_negctrl_design_matrix),
                                                 linear_terms_mat = liner_combo_vec,
                                                 type = "logit"))

    # Constructs the data to return
    #print("DEBUG gen_posctrl_data(): Creating positive control data")
    posctrl_data <- negctrl_data
    posctrl_data <- posctrl_data %>%
        mutate(!!get_outcome_colname() := postctrl_outcomes_vec)
    posctrl_data <- posctrl_data %>%
        mutate(!!get_outcomedatatype_colname_prefix() := rep("posctrl", times = nrow(posctrl_data)))

    return(posctrl_data)

} # End gen_posctrl_data()

#' From help page of sample() function that can handle passing of a single
#' number.
resample <- function(x, ...) x[sample.int(length(x), ...)]

verify_binary_outcome <- function(data = NULL, trt_var = "Z", out_var = "Y",
    expected_trt_eff = log(1.0)) {

    if (!require(dplyr)) { library(dplyr) }
    if (!require(broom)) { library(broom) }

    glm_formula_rhs <- data %>% select(-.data[[out_var]]) %>%
         colnames() %>% paste(collapse = " + ")
    glm_formula <- paste(out_var, glm_formula_rhs, sep = " ~ ")

    #print(glm_formula); browser()

    glm_fit <- glm(glm_formula, data = data, family = binomial(link = "logit"))

    glm_fit_est <- broom::tidy(glm_fit)

    trt_eff_est <- glm_fit_est %>% filter(term == trt_var) %>% pull(estimate)

    return(isTRUE(all.equal(trt_eff_est, expected_trt_eff, tolerance = 0.1)))

} # End verify_outcome()

gen_additive_pscore_formula <- function(varnames = vector("character", 0),
    trt_var = "Z", outcome_var = "Y", vars_to_remove = vector("character", 0)) {

    if (!require(vctrs)) { library(vctrs) }

    vars_to_remove_from_formula_rhs <- vec_c(trt_var, outcome_var, vars_to_remove)
    covars <- varnames[!varnames %in% stringr::coll(vars_to_remove_from_formula_rhs)]

    formula_rhs <- paste0(covars, collapse = "+")
    pscore_formula <- paste(trt_var, formula_rhs, sep = " ~ ")
    return(pscore_formula)
}

get_pscore_weights <- function(trt_data = NULL,
    pscore_formula = vector("character", 0),
    stabilised = TRUE, trim_weights_at_99pc = TRUE, trim_lower = FALSE) {

    if (!require(WeightIt)) { library(WeightIt) }

    pscores_obj <- NULL
    if (FALSE == stabilised) {
        pscores_obj <- weightit(as.formula(pscore_formula), data = trt_data, method = "ps")
    } else {
        pscores_obj <- weightit(as.formula(pscore_formula), data = trt_data,
            method = "ps", stabilize = TRUE)
    }

    if (TRUE == trim_weights_at_99pc) {
        trimmed_weights_obj <- WeightIt::trim(pscores_obj, at = 0.99, lower = trim_lower)
        return(trimmed_weights_obj$weights)
        #return(trim(pscores_obj, at = 0.99)$weights)
    } else {
        return(pscores_obj$weights)
    }
} # End get_pscore_weights()

est_trtcoef_w_sweights <- function(data = NULL, trt_var = "Z", outcome_var = "Y",
    vars_to_remove_from_trt_model = vector("character", 0),
    vars_to_remove_from_outcome_model = vector("character", 0),
    different_est_stderr_methods = FALSE) {

    if (!require(WeightIt)) { library(WeightIt) }
    if (!require(survey))   { library(survey)   }
    if (!require(glue))     { library(glue)   }

    pscore_formula <- gen_additive_pscore_formula(
        varnames = colnames(data), trt_var = trt_var, outcome_var = outcome_var,
        vars_to_remove = vars_to_remove_from_trt_model)

    print(glue(">>> DEBUG est_trtcoef_w_sweights() Treatment model `{pscore_formula}`"))
    stabilised_ipw <- get_pscore_weights(trt_data = data,
        pscore_formula = pscore_formula, stabilised = TRUE, trim_weights_at_99pc = TRUE)

    # Outcome
    varnames <- colnames(data)
    vars_to_remove_from_outcome_formula_rhs <- c(outcome_var, vars_to_remove_from_outcome_model)
    outcome_model_covars <- varnames[!varnames %in% stringr::coll(vars_to_remove_from_outcome_formula_rhs)]

    outcome_formula_rhs <- paste0(outcome_model_covars, collapse = "+")
    outcome_model_formula <- paste(outcome_var, outcome_formula_rhs, sep = " ~ ")
    print(glue("> DEBUG: est_trtcoef_w_sweights() Outcome model `{outcome_model_formula}`"))

    if (different_est_stderr_methods) {
        est_trt_effect <- create_trt_eff_estimator_with_formula(outcome_model_formula)
        logit_fit_res <- est_trt_effect(data = data, weights = stabilised_ipw) %>%
            broom::tidy()
    }

    # Sandwich estimator CI
    logit_fit_survey <- svyglm(as.formula(outcome_model_formula),
        family = quasibinomial(link = logit),
        design = svydesign(id = ~ 1,                   # No cluster
            weights = stabilised_ipw, data = data))    # IPW weights
    logit_fit_survey_res <- broom::tidy(logit_fit_survey)

    if (!different_est_stderr_methods) {
        logit_fit_res <- logit_fit_survey_res
    }

    if (different_est_stderr_methods) {
        # Obtains SE from `sandwich estimator`
        logit_fit_res <- logit_fit_res %>% mutate(
            std.error = logit_fit_survey %>% tidy() %>% pull(std.error))
    }

    return(logit_fit_res %>% select(c(term, estimate, std.error)) )
} # End `est_trtcoef_w_sweights()`

convert_negctrl_est_to_coef_tbl <- function(negctrl_est_tbl = NULL,
    num_negctrl = 0L, trt_covarname = "Z", confounder_est_prefix = "X",
    new_coef_prefix = "c") {

    if (!require(dplyr))     { library(dplyr)   }
    if (!require(stringr))   { library(stringr)   }

    # Assumes `negctrl_est_tbl` is in the form
    # outcome_id,  est
    #          2,  list<GLM est>

    postctrl_coef_mat <- negctrl_est_tbl %>% apply(MARGIN = 1, FUN = function(row) {

        outcome_id <- row[[1]]
        negctrl_est <- row[[4]]

        negctrl_est_est <- negctrl_est %>% select(term, estimate)

        # Replaces the column names
        negctrl_colnames <- negctrl_est_est %>% pull(term)

        coef_colnames <- negctrl_colnames %>% stringr::str_replace(
            pattern = "\\(Intercept\\)",
            replace = paste0(new_coef_prefix, 0))

        coef_colnames <- coef_colnames %>% stringr::str_replace(
            pattern = paste0(trt_covarname),
            replace = paste0(new_coef_prefix, trt_covarname))

        coef_colnames <- coef_colnames %>% stringr::str_replace(
            pattern = paste0("^", confounder_est_prefix),
            replace = new_coef_prefix)

        # Adds outcome ID covarname
        coef_colnames <- c("outcome_id", coef_colnames)

        coef_to_rt <- negctrl_est_est %>% pull(estimate)
        # Adds outcome ID
        coef_to_rt <- c(outcome_id, coef_to_rt)

        names(coef_to_rt) <- coef_colnames
        return(coef_to_rt)

    }) # End of apply()

    postctrl_coef_tbl <- postctrl_coef_mat %>% t() %>% as_tibble()
    postctrl_coef_tbl[['outcome_id']] <- as.integer(postctrl_coef_tbl[['outcome_id']])

    return(postctrl_coef_tbl)

} # End convert_negctrl_est_to_coef_tbl()

#' Estimates regression coefficients of positive controls
est_coef_posctrl <- function(outcome_id = NA_integer_, negctrl_id = NA_integer_,
                             postctrl_data = NULL,
                             outcome_varname = "Y",
                             trt_varname = "Z", use_ipw = TRUE,
                             vars_to_remove_from_trt_model = "",
                             vars_to_remove_from_outcome_model = "") {

    if (!require(dplyr))   { library(dplyr)   }
    if (!require(stringr)) { library(stringr) }
    if (!require(broom))   { library(broom)   }

    if (is.null(postctrl_data)) { error("!! ERROR: est_coef_posctrl() Positive control data is NULL") }

    print(glue("> DEBUG: est_coef_posctrl() Fitting IPW GLM to +POSITIVE+ control ",
               # "true effect size {signif(trt_eff, 3)}, generated using ",
               "neg control {negctrl_id} (outcome num {outcome_id})"))

    # Manual GLM
    # covars_in_data <- colnames(postctrl_data)
    # covar_terms <- covars_in_data[covars_in_data != outcome_varname]
    # glm_formula_rhs <- paste0(covar_terms, collapse = " + ")
    #
    # glm_formula <- paste(outcome_varname, glm_formula_rhs, sep = "~")
    #
    # logit_fit <- glm(formula = as.formula(glm_formula),
    #                  family = binomial(link = "logit"),
    #                  data = postctrl_data)

    posctrl_coef_est_tbl <- est_trtcoef_w_sweights(postctrl_data,
                           trt_var = trt_varname,
                           outcome_var = outcome_varname,
                           vars_to_remove_from_trt_model = vctrs::vec_c(trt_varname, vars_to_remove_from_trt_model),
                           vars_to_remove_from_outcome_model = vctrs::vec_c(outcome_varname, vars_to_remove_from_outcome_model))

    # Makes sure treatment estimate is close to what was expected.

    return(list(outcome_id = outcome_id, negctrl_id = negctrl_id,
                est = list(posctrl_coef_est_tbl)) )

} # End est_coef_posctrl()


#' Collates the estimates from negative and positive controls to create data
#' structure for empirical calibration.
#'
#' @param negctrl_est_tbl \code{data.frame} with columns in the order of
#'     \code{outcome_id, negctrl_id, true_trt_eff, est}, with \code{est}
#'     a \code{list data.frame} with regression estimates in tidy format.
#' @param posctrl_est_tbl \code{data.frame} with columns in the order of
#'     \code{outcome_id, negctrl_id, true_trt_eff, est}, with \code{est}
#'     a \code{list data.frame} with regression estimates in tidy format.
#' @param trt_varname Variable name of treatment
collate_est_from_ctrls <- function(negctrl_est_tbl = NULL,
                                   posctrl_est_tbl = NULL,
                                   trt_varname = "Z",
                                   posctrl_has_adj_true_eff = FALSE,
                                   ctrl_tbl_colnames = list(
                                       outcome_id_colname = 'outcome_id',
                                       negctrl_id_colname = 'negctrl_id',
                                       est_colname = 'est',
                                       true_trt_eff_colname = 'true_trt_eff',
                                       adj_true_trt_eff_colname = 'adj_true_trt_eff')) {

    if (!require(dplyr))   { library(dplyr) }

    if (is.null(negctrl_est_tbl) | is.null(posctrl_est_tbl)) {
        error("!! ERROR: collate_est_from_ctrls() Estimates data structure is NULL")
    }
    extract_est_row <- function(row, trt_varname = "Z", has_adj_true_eff = FALSE) {

        true_eff <- row[[ ctrl_tbl_colnames[['true_trt_eff_colname']] ]]
        if (has_adj_true_eff) {
            adj_true_eff_vec <- row[[ ctrl_tbl_colnames[['adj_true_trt_eff_colname']] ]]
        } else {
            adj_true_eff_vec <- true_eff
        }

        actual_est <- row[[ ctrl_tbl_colnames[['est_colname']] ]] %>%
            filter(term == trt_varname) %>% pull(estimate)
        actual_est_stderr <- row[[ ctrl_tbl_colnames[['est_colname']] ]] %>%
            filter(term == trt_varname) %>% pull(std.error)

        new_row <- c(outcome_id = row[[ ctrl_tbl_colnames[['outcome_id_colname']] ]],
                 negctrl_id = row[[ ctrl_tbl_colnames[['negctrl_id_colname']] ]],
                 true_eff = true_eff,
                 actual_eff_est = actual_est,
                 actual_eff_est_stderr = actual_est_stderr,
                 adj_true_eff = adj_true_eff_vec)

        return (new_row)
    }

    negctrl_est <- negctrl_est_tbl %>% apply(MARGIN = 1, FUN = extract_est_row,
                                             has_adj_true_eff = FALSE)

    posctrl_est <- posctrl_est_tbl %>% apply(MARGIN = 1, FUN = extract_est_row,
                                             has_adj_true_eff = posctrl_has_adj_true_eff)

    est_to_rt <- cbind(negctrl_est, posctrl_est) %>% t() %>% as_tibble()

    return(est_to_rt)

} # End collate_est_from_ctrls()

#' Extract results from simulation
extract_tbl_from_sim_result <- function(sim_res = NULL, name = "") {

    if (!require(purrr)) { library(purrr) }

    res_tbl <- sim_res %>% purrr::map_dfr(.f = ~{
        result <- .x
        calib_result <- result[[name]]
        return(list(sim_num = 0L, res = list(calib_result)))
    })

    res_tbl <- res_tbl %>% mutate(sim_num = seq.int(from = 1L, to = length(sim_res)))

} # End extract_from_sim_result()

#' Merges calibrated and uncalibrated controls data.
merge_calib_uncalib_ctrl_data <- function(uncalibrated_ctrl_data,
                                          calibrated_ctrl_data,
                                          cali_sys_err_model_type = "nomod") {

    calibrated_estimates <- calibrated_ctrl_data %>%
        mutate(sim_id = uncalibrated_ctrl_data[['sim_id']]) %>%
        mutate(outcome_id = uncalibrated_ctrl_data[['outcome_id']]) %>%
        mutate(negctrl_id = uncalibrated_ctrl_data[['negctrl_id']]) %>%
        mutate(true_eff  = uncalibrated_ctrl_data[['true_eff']]) %>%
        relocate(true_eff) %>% relocate(negctrl_id) %>% relocate(outcome_id) %>% relocate(sim_id) %>%
        mutate(restype = rep("calibrated", times = nrow(calibrated_ctrl_data))) %>%
        mutate(calitype = rep(NA_character_, times = nrow(calibrated_ctrl_data)))

    ctrl_data_to_calibrate <- uncalibrated_ctrl_data %>%
        mutate(est_L95 = actual_eff_est  - 1.95 * (actual_eff_est_stderr )) %>%
        relocate(est_L95, .after = actual_eff_est) %>%
        mutate(est_U95 = actual_eff_est  + 1.95 * (actual_eff_est_stderr )) %>%
        relocate(est_U95, .after = est_L95) %>%
        mutate(restype = rep("uncalibrated", times = nrow(uncalibrated_ctrl_data))) %>%
        mutate(calitype = rep(cali_sys_err_model_type, times = nrow(calibrated_ctrl_data)))

    common_cols <- intersect(names(calibrated_estimates), names(ctrl_data_to_calibrate))

    plot_data <- rbind(ctrl_data_to_calibrate[common_cols], calibrated_estimates[common_cols])
} # End merge_calib_uncalib_ctrl_data()

#' Filters simulations -without- NA in a column.
#'
#' @param data \code{data.frame} to check
#' @param na_col Name of column to check for NA value.
#' @param group_var Column to group estimates.
filter_sims_without_na <- function(data,
    na_col = "actual_eff_est_stderr", group_var = "sim_num") {

    unique_groups_with_na <- data %>% filter(is.na(.data[[na_col]])) %>%
            pull(.data[[group_var]]) %>% unique()

    num_unique_groups_with_na <- unique_groups_with_na %>% length()

    data_without_na <- data %>%
        filter(!(.data[[group_var]] %in% unique_groups_with_na))

    return(list(data_without_na = data_without_na,
        na_groups = unique_groups_with_na))

} # End filter_sims_without_na()

#' Calculates standard deviation (sigma) of a Gaussian distribution given
#' the mean and the mean deviation of samples.
#' See https://mathworld.wolfram.com/MeanDeviation.html
calc_gaussian_sd_from_mean_deviation <- function(mu = 0.0, mean_deviation = 0.2) {
	var_target <- (pi / 2) * (mean_deviation)^2
	sd_target <- sqrt(var_target)
	return(sd_target)
} # End calc_gaussian_sd_from_mean_deviation()
