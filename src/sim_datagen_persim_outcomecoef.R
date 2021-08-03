gen_outcome_model_coef_persim <- function(
    coef_in_primaryoutcome = TRUE, coef_in_negctrls = TRUE,
    coef_range = c(min_or = 0.5, max_or = 2.0),
    num_outcomes = 2L,
    primary_outcome_num = 1L,
    coef_noeff_val = NA_real_) {

    if (!require(glue)) { library(glue) }

    if (num_outcomes < 2) {
        warning(glue("!! gen_outcome_model_coef_persim(): Required more than two outcomes."))
    }

    beta_vec <- rep(coef_noeff_val, times = num_outcomes)
    if (isTRUE(coef_in_primaryoutcome)) {
        # Outcome of interest
        beta_vec[[primary_outcome_num]] <-
            log(runif(1, min = coef_range[['min_or']], max = coef_range[['max_or']]))
    }

    if (isTRUE(coef_in_negctrls)) {
        num_negctrls <- num_outcomes - 1L
        # Negative control outcomes depend on quad term.
        beta_vec[-primary_outcome_num] <-
            log(runif(num_negctrls, min = coef_range[['min_or']], max = coef_range[['max_or']]))
    }

    return(beta_vec)
} # End gen_outcome_model_coef_persim()

gen_same_outcome_model_coef_persim <- function(
    coef_range = c(min_or = 0.5, max_or = 2.0),
    num_outcomes = 2L,
    coef_noeff_val = NA_real_) {

    beta_vec <- rep(coef_noeff_val, times = num_outcomes)

    beta_all_outcomes <- log( runif(1,
                                   min = coef_range[['min_or']],
                                   max = coef_range[['max_or']]))

    beta_vec_final <- set_coef_vals(
        initial_coef_vec = beta_vec,
        prop_negctrl_to_set = 1.0,
        val = beta_all_outcomes,
        coef_noeff_val = coef_noeff_val)

    return(beta_vec_final)
}
