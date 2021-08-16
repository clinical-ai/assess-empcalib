suppressPackageStartupMessages({
    library(tidyr)
    library(dplyr)
    library(tibble)
    library(glue)
    library(here)
    library(purrr)
    library(vctrs)
    library(broom)
    library(ggplot2)
    library(readr)

    library(assertthat)

    library(doParallel)

    library(cobalt)    # For propensity balance plotting
    library(WeightIt)  # For obtaining IPW
    library(survey)    # For obtaining CI from `sandwich` estimator, aka robust estimator

    library(EmpiricalCalibration)

    # For R-3.6.3
    if (as.numeric(unlist(getRversion()))[[1]] < 4) {
        library(backports)
    }
})

source(here::here("src", "fn_prog_helpers.R"))
source(here::here("src", "sim_datagen.R"))
source(here::here("src", "sim_datagen_persim_outcomecoef.R"))
source(here::here("src", "sim_binaryoutcome_utils.R"))
source(here::here("src", "plot_funnel_est_coverage.R"))
source(here::here("src", "plot_bias_boxplot-v2arch.R"))
source(here::here("src", "plot_forestplots-v2arch.R"))
source(here::here("src", "plot_ctrl_empcalib_forestplot.R"))
source(here::here("src", "calib_ctrl.R"))
source(here::here("src", "calib_posctrls.R"))
source(here::here("src", "read_configfile.R"))

# Reads in setting file if supplied ----------------

args <- commandArgs(trailingOnly = TRUE)
settings_file <- args[1]  # NA if not supplied.
if (is.na(settings_file)) {
    settings_file <- here::here("data", "config", "unmeasuredconf",
        "simconfig-unmeasuredconf-allrand-5negctrl.csv")
    print(glue("> DEBUG: Setting file not supplied. Reading DEFAULT from `{settings_file}`"))
}
print(glue("> DEBUG: Reading settings from `{settings_file}`"))
get_config <- load_config(csvconfigfile = settings_file,
    configcols = list(setting_colname = "flag", value_colname = "value"))

n_obs <- get_config("n_obs_persim", default_val = 10000) %>% as.integer()

# n_sim ----
n_sim <- ifelse(!is.na(settings_file),
                yes = get_config("n_sims", default_val = 0L), no = 10L) %>% as.integer()
parallelise_sim_threshold <- 20L

# Interrupt exec flow at debug points, reduces the number of simulations
# below `parallelise_sim_threshold`
debug_intrpt <- ifelse(!is.na(settings_file),
                        yes = get_config("debug_intrpt", default_val = "FALSE") %>% as.logical(),
                        no = FALSE)
if (debug_intrpt) { n_sim <- 1L }

# Number of negative controls
num_negctrls <- ifelse(!is.na(settings_file),
                       yes = get_config("n_negctrls", default_val = 5L),
                       no = 5L) %>% as.integer()
num_notpriout_mconf <- ifelse(!is.na(settings_file),
                              yes = get_config("num_not_trt_pri_outcome_mconfounder", default_val = 5L),
                              no = 5L) %>% as.integer()
num_outcomes <- num_negctrls + 1L # Plus one for outcome of interest

# Number of processes to launch
nprocesses <- ifelse(!is.na(settings_file),
                     yes = get_config("nprocess", default_val = 4L) %>% as.integer(),
                     no = 4L)

# Press ENTER or a key to manually move to next plot.
manual_forward_plots <- FALSE

# Use the same treatment model coef. across all simulations
all_sim_use_trt_model_template <- FALSE

# Persists simulation data
persist_sim_out_data <- TRUE
persist_plots <- TRUE

# Simulation number to persist plots.
sim_num_to_plot <- as.integer(n_sim / 2)

# If RANDOM coefficients
seed_rand <- ifelse(!is.na(settings_file),
                    yes = get_config("seed_outcome_coef_rand", default_val = "FALSE") %>% as.logical(),
                    no = FALSE)
if (seed_rand) { set.seed(12890) }

# Value of systematic error model
sys_err_model_legacy_flag <- FALSE

# Name of simulation. Used as part of prefixes of filenames. # Name of simulation ----
sim_shortdesc <- ifelse(!is.na(settings_file),
                        yes = get_config("sim_shortdesc", default_val = "TEST"),
                        no = "TEST") %>% as.character()
sim_name_part <- glue("{sim_shortdesc}-n{n_obs}-s{n_sim}-nctrl_{num_negctrls}-n_mconf_{num_notpriout_mconf}")

no_eff_loghr <- log(1.0)
weak_trt_eff_loghr <- log(1.25)
mid_trt_eff_loghr <- log(1.5)
strong_trt_eff_loghr <- log(1.75)
vstrong_trt_eff_loghr <- log(2.0)
extreme_strong_trt_eff_loghr <- log(4.0)

# Positive control trt effects ----
posctrl_log_or_targets <- log( c(2.0, 4.0, 6.0))

# Calibrates negative controls similar to Suchard, Schuemie, Lancet 2019
calibrate_negctrls <- TRUE
calibrate_posctrls <- TRUE

##
# Simulation types ------------------------------------------------------------
##

set_binary_status <- function(status) { return(function() { status }) }

set_value <- function(prop) {
    prop_negctrl_same_coef <- prop
    return(
        function() { return(prop_negctrl_same_coef) }
    )
}

## Flags for unmeasured confounding scenario ----

### Extra confounder, U, in study only?
gen_extra_confounder_in_study <- get_config(
    "gen_extra_confounder_in_study", default_val = FALSE) %>% as.logical()
    
### Extra confounder, U, in negative controls?
gen_extra_confounder_in_negctrl <- get_config(
    "gen_extra_confounder_in_negctrl", default_val = FALSE) %>% as.logical()

### U term becomes unmeasured confounder in outcome of interest
ignore_extra_confounder_in_study <- get_config(
    "ignore_extra_confounder_in_study", default_val = FALSE) %>% as.logical()

### U term becomes unmeasured confounder in negative controls
ignore_extra_confounder_in_negctrl <- get_config(
    "ignore_extra_confounder_in_negctrl", default_val = FALSE) %>% as.logical()

# Proportion of outcomes affected by this potential unmeasured confounder `U`
prop_beta_u_with_effect <- 1.0

sim_name_part <- paste(sim_name_part,
                       ifelse(gen_extra_confounder_in_negctrl,
                              yes = "nctrlcap_y", no = "nctrlcap_n"),
                       sep = "-")

negctrl_extra_confounder_coef_relative_outcome_interest <- set_binary_status(FALSE)

static_extra_confounder_coef <- set_binary_status(FALSE)
#sim_name_part <- paste(sim_name_part,
#                       ifelse(gen_extra_confounder_in_negctrl,
#                              yes = "static_u_y", no = "static_u_n-"),
#                       sep = '-')

## Flags for model mis-specification of quadratic term scenario ----

### main flag indicating the scenario
quad_term_sim_scenario     <- ifelse(!is.na(settings_file),
                                     yes = get_config("quad_term_sim_scenario", default_val = FALSE) %>% as.logical(),
                                     no = FALSE)
x1sqr_in_trt_model         <- NA_integer_
gen_study_with_x1sqr       <- NA_integer_
gen_negctrls_with_x1sqr    <- NA_integer_
# Flags for estimation
ignore_x1sqr_est_trt       <- NA_integer_
ignore_x1sqr_study_model   <- NA_integer_
ignore_x1sqr_negctrl_model <- NA_integer_

set_status_for_alltrue_deps <- function(dep_cond_bool_vec = c()) {
    return(function(status) {
        status & ( sum(TRUE == (dep_cond_bool_vec)) == length(dep_cond_bool_vec))
    })
}

if (isTRUE(quad_term_sim_scenario)) {
    # Presence/absence of quad term when generating treatment
    x1sqr_in_trt_model         <- TRUE # Always true

    # Presence/absence of quad term in generating outcomes
    gen_study_with_x1sqr       <- TRUE
    gen_negctrls_with_x1sqr    <- ifelse(!is.na(settings_file),
                                     yes = get_config("gen_negctrls_with_x1sqr", default_val = FALSE) %>% as.logical(),
                                     no = FALSE)

    set_logical_primary_outcome_with_x1sqr <- set_status_for_alltrue_deps(
        c(quad_term_sim_scenario, gen_study_with_x1sqr))
    set_logical_negctrl_outcome_with_x1sqr <- set_status_for_alltrue_deps(
        c(quad_term_sim_scenario, gen_negctrls_with_x1sqr))

    # Presence/absence of quad term in OUTCOME model (estimation). Ignoring quad term introduces bias.
    ignore_x1sqr_est_trt       <- ifelse(FALSE == x1sqr_in_trt_model, yes = FALSE,
                                         no = TRUE)
    ignore_x1sqr_study_model   <- set_logical_primary_outcome_with_x1sqr(TRUE)
    ignore_x1sqr_negctrl_model <- set_logical_primary_outcome_with_x1sqr(TRUE)
}

## Flags for scenario: Interaction of two confounders X1X2 -----
## Requires at least one negative control and two confounders

intend_intr_twoconfounders_sim_scenario <- ifelse(!is.na(settings_file),
    yes = get_config("intend_intr_twoconfounders_sim_scenario", default_val = FALSE) %>% as.logical(),
    no = FALSE)
# Conditional on number measured confounders
intr_twoconfounders_sim_scenario <- ((num_notpriout_mconf >= 2) & intend_intr_twoconfounders_sim_scenario)

# Initialise settings
gen_study_with_x1x2 <- NA
gen_negctrls_with_x1x2 <- NA
ignore_x1x2_in_obsest <- NA

if (intr_twoconfounders_sim_scenario) {
    # Presence/absence of X1X2 term when generating data
    gen_study_with_x1x2      <- TRUE
    gen_negctrls_with_x1x2   <- ifelse(!is.na(settings_file),
                                     yes = get_config("gen_negctrls_with_x1x2", default_val = FALSE) %>% as.logical(),
                                     no = FALSE)
    # Presence/absence of X1X2 term in observed data
    ignore_x1x2_in_obsest    <- TRUE
}

## Flags for model mis-specification: Interaction of confounder and treatment -----

intr_conftrt_sim_scenario     <- FALSE    # main flag indicating the scenario

# Initial values of settings
ignore_intr_conftrt_in_obsest <- NA

if (isTRUE(intr_conftrt_sim_scenario)) { # settings relevant to this scenario
    # Presence/absence of X1*Z term in observed data
    ignore_intr_conftrt_in_obsest <- TRUE

} # End TRUE == intr_conftrt_sim_scenario

## Flags for lack of positivity -----

# main flag indicating the scenario
lackpositivity_sim_scenario <- ifelse(
    !is.na(settings_file),
    yes = get_config("lackpositivity_sim_scenario", default_val = FALSE) %>% as.logical(),
    no = FALSE)

trt_assign_var_threshold <- 0.65

lackpositivity_coef_inflation_fct <- ifelse(
    !is.na(settings_file),
    yes = get_config("lackpositivity_coef_inflation_fct", default_val = 1.20) %>% as.numeric(),
    no = 1.20)
if (lackpositivity_coef_inflation_fct < 1.0) {
    warning(glue(
        "!! WARN: Inflation factor `lackpositivity_coef_inflation_fct` ",
        "{signif(lackpositivity_coef_inflation_fct, 4)} is less than 1.0; will be DEFLATING coef."))
}

if (isTRUE(lackpositivity_sim_scenario)) {
    gen_negctrls_with_trtpos <- ifelse(
        !is.na(settings_file),
        yes = get_config("gen_negctrls_with_trtpos", default_val = FALSE) %>% as.logical(),
        no = FALSE)
}

## Flags for measurement error of a confounder ----

measurement_err_confounder_scenario <- ifelse(
    !is.na(settings_file),
    yes = get_config("merror_confounder_scenario", default_val = FALSE) %>% as.logical(),
    no = FALSE)

gen_negctrls_with_merror_covar <- ifelse(
    !is.na(settings_file),
    yes = get_config("gen_negctrls_with_merror_covar", default_val = FALSE) %>% as.logical(),
    no = FALSE)

rand_merror_gaussian_params <- ifelse(
    !is.na(settings_file),
    yes = get_config("rand_merror_gaussian_params", default_val = FALSE) %>% as.logical(),
    no = FALSE)

merror_coef_oddsratio_inflation_fct <- NA_real_
rand_merror_gaussian_mean_min <- NA_real_
rand_merror_gaussian_mean_max <- NA_real_
rand_merror_gaussian_mean_deviation_min_pct <- NA_real_
rand_merror_gaussian_mean_deviation_max_pct <- NA_real_
merror_gaussian_mean <- NA_real_
merror_gaussian_sd <- NA_real_

if (isTRUE(measurement_err_confounder_scenario)) {
    merror_coef_oddsratio_inflation_fct <- ifelse(
            !is.na(settings_file),
            yes = get_config("merror_coef_oddsratio_inflation_fct", default_val = FALSE) %>% as.numeric(),
            no = 1.20) %>%
        if_cond_do(condition = (lackpositivity_coef_inflation_fct < 1.0),
            fun = { function(inflate_scale) {
                warning(glue(
                    "!! WARN: Measurement error inflation factor `merror_coef_oddsratio_inflation_fct` ",
                    "{signif(merror_coef_oddsratio_inflation_fct, 4)} is less than 1.0; will be set to 1.0"))
                return(1.0)
            } })
}

if (measurement_err_confounder_scenario & rand_merror_gaussian_params) {
    # Randomise mean and standard error of error
    rand_merror_gaussian_mean_min <- ifelse(
        !is.na(settings_file),
        yes = get_config("rand_merror_gaussian_mean_min", default_val = 1.0) %>% as.numeric(),
        no = 1.0)
    rand_merror_gaussian_mean_max <- ifelse(
        !is.na(settings_file),
        yes = get_config("rand_merror_gaussian_mean_max", default_val = 2.0) %>% as.numeric(),
        no = 2.0)
    rand_merror_gaussian_mean_deviation_min_pct <- ifelse(
        !is.na(settings_file),
        yes = get_config("rand_merror_gaussian_mean_deviation_min_pct", default_val = 10.0) %>% as.numeric(),
        no = 10.0)
    rand_merror_gaussian_mean_deviation_max_pct <- ifelse(
        !is.na(settings_file),
        yes = get_config("rand_merror_gaussian_mean_deviation_max_pct", default_val = 20.0) %>% as.numeric(),
        no = 20.0)
} # End if measurement_err_confounder_scenario & isTRUE(rand_merror_gaussian_params)

if (measurement_err_confounder_scenario & !rand_merror_gaussian_params) {
    # Obtain non-randomised mean and standard error of error
    merror_gaussian_mean <- ifelse(
        !is.na(settings_file),
        yes = get_config("merror_gaussian_mean", default_val = 0.0) %>% as.numeric(),
        no = 0.0)
    merror_gaussian_mean <- merror_gaussian_mean %>%
        if_cond_do(condition = merror_gaussian_mean < 0.0,
            fun = { function(inflate_scale) {
                warning(glue(
                    "!! WARN: Measurement Gaussian mean `merror_gaussian_mean`
                     {signif(merror_gaussian_mean, 4)} is less than 0.0; will be set to 0.0"))
                return(0.0)
            } })

    merror_gaussian_sd <- ifelse(
        !is.na(settings_file),
        yes = get_config("merror_gaussian_sd", default_val = 1.0) %>% as.numeric(),
        no = 1.0)
    merror_gaussian_sd <- merror_gaussian_sd %>%
        if_cond_do(condition = (merror_gaussian_sd < 0.0),
            fun = { function(inflate_scale) {
                warning(glue(
                    "!! WARN: Measurement Gaussian SD `merror_gaussian_sd`
                     {signif(merror_gaussian_sd, 4)} is less than 0.0; will be set to 1.0"))
                return(1.0)
            } })
} # End if measurement_err_confounder_scenario & !rand_merror_gaussian_params

## Other settings ----

# beta_u_s , effect size for negative controls, is relative to beta_u from outcome of interest.


static_extra_confounder_coef_all_sim_iters <- set_binary_status(FALSE)
b_u_all_sim_iters <- NA_real_  # in log-odds scale.
if (static_extra_confounder_coef_all_sim_iters()) {
    b_u_all_sim_iters <- log(4.0)
}

sim_name <- sim_name_part # The latest content becomes the simulation name.

##
# Generate the confounder specification ----------------------------------------
##
print(glue("> DEBUG: Total number of neg ctrl outcomes={num_negctrls}, ",
           " num independent confounders={num_notpriout_mconf}"))
confounder_spec <- create_confounder_spec(n = num_notpriout_mconf + 1L,
                                          datagen_fn = "rnorm",
                                          confounder_name_prefix = get_confounder_colname_prefix())

potential_unmeasured_confounder_varname <- paste0(get_extra_confounder_colname(), 1)

if (gen_extra_confounder_in_study | gen_extra_confounder_in_negctrl) {
    # Generate values for potential unmeasured confounder
    print("> DEBUG: Adding spec. to generate potential unmeasured confounder U.")
    confounder_spec <- confounder_spec %>% add_one_confounder_spec(
        datagen_fn = "rnorm", genfn_args = list(mean = 0),
        covar_name = potential_unmeasured_confounder_varname)
}

if (lackpositivity_sim_scenario) {
    # Specifies additional covariate to generate to determine lack of positivity
    print("> DEBUG: Adding spec. to generate extra confounder for lack of positivity scenario")
    confounder_spec <- confounder_spec %>% add_one_confounder_spec(
        datagen_fn = "rnorm", genfn_args = list(mean = 0),
        covar_name = get_trtpos_deltavar_colname())
}

if (measurement_err_confounder_scenario) {
    # Specifies additional covariate to generate that will eventually have
    # measurement error
    print("> DEBUG: Adding spec. to generate extra confounder for measurement error scenario")
    confounder_spec <- confounder_spec %>% add_one_confounder_spec(
        datagen_fn = "rnorm", genfn_args = list(mean = 0),
        covar_name = get_merror_confounder_colname())
} # End add to confounder_spec for measurement_err_confounder_scenario

if (debug_intrpt) {
    print("> DEBUG: Spec. for independent confounders generated. `confounder_spec`")
    print(confounder_spec) ; browser()
}

##
# Build templates of regression coefficients ----------------------------------
##

###
## DataGen, Treatment model: template for coefficients in treatment model -------------------------
###

print(glue("> DEBUG: Generating regression coefficients TEMPLATE for ",
    "TREATMENT MODEL, with {num_negctrls} negative controls and ",
    "{num_notpriout_mconf} independent fully measured confounders."))

alpha_coef_range <- get_coef_range_or()

trt_model_coefs_template_vec <- confounder_spec %>% gen_trt_model_regression_coefs(
    coef_range = alpha_coef_range)

if (all_sim_use_trt_model_template) {
    print("> DEBUG: Using same TREATMENT MODEL for ALL simulations")
}

### DataGen, Treatment model template: Quad term X1Sqr -----

alpha_x1sqr_range <- get_coef_range_or()
# If required, models the primary confounder X1 to treatment association
# as quadratic, by adding coefficient for squared term X1Sqr.
# Updates the coefficients in treatment model for quadratic term.
#
# Need to to this because quad term depends on first covariate, hence not included
# in initial `confounder_spec`
if (quad_term_sim_scenario & x1sqr_in_trt_model) {
    print("> DEBUG: Adding coefficient for quad term in TREATMENT model")

    old_names <- names(trt_model_coefs_template_vec)
    trt_model_coefs_template_vec <- c(trt_model_coefs_template_vec,
                                      log(runif(1,
                                                min = alpha_x1sqr_range[['min_or']],
                                                max = alpha_x1sqr_range[['max_or']])) )
    names(trt_model_coefs_template_vec) <- c(old_names, get_quadterm_colname())
}

### DataGen, Treatment model template: X1X2 interaction -----

if (intr_twoconfounders_sim_scenario) {
    print("> DEBUG: Adding coefficient of X1X2 term in TREATMENT model TEMPLATE")

    alpha_x1x2_range <- get_coef_range_or()

    old_names <- names(trt_model_coefs_template_vec)
    trt_model_coefs_template_vec <- c(trt_model_coefs_template_vec,
                                      log(runif(1,
                                                min = alpha_x1x2_range[['min_or']],
                                                max = alpha_x1x2_range[['max_or']])) )
    names(trt_model_coefs_template_vec) <- c(old_names, get_x1x2intr_colname())
}

### DataGen, Treatment model template: Lack of positivity case -----
if (lackpositivity_sim_scenario) {
    print("> DEBUG: Inflating coefficient of covar to adjust positivity in treatment model")

    beta_xpos_deltavar_range <- get_coef_range_or()

    # Ensures this coefficient is always highest
    xtrt_log_or <- log(beta_xpos_deltavar_range[['max_or']] * lackpositivity_coef_inflation_fct)
    trt_model_coefs_template_vec[[get_trtpos_deltavar_colname()]] <- xtrt_log_or
}

###
## DateGen: Build outcome model coef TEMPLATE ------------------------------------------
###

# True treatment effect in outcome of interest.
trt_eff_range <- get_coef_range_or()
trt_eff <- log( runif(1,
                min = trt_eff_range[['min_or']], max = trt_eff_range[['max_or']]))

# Use the same outcome model template, for measured confounders,
# across ALL simulations.
use_outcome_model_coef_template_all_sims <- FALSE

# Coefficients of measured confounders
measured_conf_coef_range <- get_coef_range_or()

# All confounders effect all outcomes with the *same magnitude*
# (`confounder_to_all_outcomes_all_sims_eff`) across all simulations. aka Strategy one.
outcome_model_coef_confounder_same_eff_all_outcomes_all_sims <- FALSE
confounder_to_all_outcomes_all_sims_eff <- log(3.5)

print(glue("> DEBUG: Generating OUTCOME MODEL regression coefficients TEMPLATE for MEASURED confounders ",
    " with {num_negctrls} negative controls, and ",
    "{num_notpriout_mconf} independent fully measured confounders."))

outcome_model_coef_prefix <- "b"

# Initial template for outcome model, containing all NAs.
outcome_model_coef_tbl_template <- gen_outcome_model_coef_template(
    total_num_negctrls = num_negctrls,
    num_not_pri_outcome_confounders = num_notpriout_mconf,
    trt_eff_logor_outcomeinterest = NA_real_,
    trt_eff_range_or = trt_eff_range,
    coef_prefix = "b", trt_valname = get_trt_colname())


if (use_outcome_model_coef_template_all_sims) {
    print("> DEBUG: Using the same outcome model template for ALL simulations")

    if (outcome_model_coef_confounder_same_eff_all_outcomes_all_sims) {
        outcome_model_coef_tbl_template <- outcome_model_coef_tbl_template %>%
            set_nontrt_varcoefs_sameval(coef_prefix = outcome_model_coef_prefix,
                                        coef_val = confounder_to_all_outcomes_all_sims_eff)
    } else {
        outcome_model_coef_tbl_template <- outcome_model_coef_tbl_template %>%
            set_nontrt_varcoefs(coef_prefix = outcome_model_coef_prefix,
                                coef_val_range = measured_conf_coef_range)

    }
} # End `TRUE == use_outcome_model_coef_template_all_sims`

# Init parameters
same_all_measured_eff_to_all_outcomes_all_sim <- NA
same_measured_eff_to_all_outcomes_per_sim_iter <- NA
# Each simulation iteration generates its own outcome model
if (FALSE == use_outcome_model_coef_template_all_sims) {
    # In each simulation,
    # coef. of ALL measured confounders, for ALL outcomes, are the same.
    same_all_measured_eff_to_all_outcomes_all_sim <- get_config(
     "same_all_measured_eff_to_all_outcomes_all_sim", default_val = FALSE) %>% as.logical()

    # In each simulation,
    # each measured confounders have the same effect size across all outcomes,
    # BUT the coef. differ across confounders
    # (unlike `same_all_measured_eff_to_all_outcomes_all_sim`)
    same_measured_eff_to_all_outcomes_per_sim_iter <- get_config(
        "same_measured_eff_to_all_outcomes_per_sim_iter", default_val = FALSE) %>% as.logical()
}

if (debug_intrpt) {
    print("> DEBUG: DataGen Outcome model template `outcome_model_coef_tbl_template`")
    print(outcome_model_coef_tbl_template) ; browser()
}

####
### DataGen, Outcome model: Outcome model coef template: Potential Unmeasured Confounder U1 --------
####

use_same_extra_term_coef_vec_all_sim <- ifelse(
    !is.na(settings_file),
    yes = get_config("use_same_extra_term_coef_vec_all_sim", default_val = FALSE) %>% as.logical(),
    no = FALSE)

# Equi-confounding (same effect of potential unmeasured confounder affects ALL outcomes)
# flag for potential unmeasured confounder U1. Applicable WITHIN and ACROSS simulations
equi_confounding <- FALSE
beta_u_for_all_outcomes <- log(2.0)

# Over-ride equi-confounding and its coefficient, when the outcome coefficient
# template setting to make *all* coefficients the same across *all* simulations
if (TRUE == outcome_model_coef_confounder_same_eff_all_outcomes_all_sims) {
    print(glue("> DEBUG: `outcome_model_coef_confounder_same_eff_all_outcomes_all_sims` ",
               "{outcome_model_coef_confounder_same_eff_all_outcomes_all_sims}. ",
               "Overriding equi-confounding setting due to outcome model template."))
    equi_confounding <- TRUE
    beta_u_for_all_outcomes <- confounder_to_all_outcomes_all_sims_eff
}

beta_u_range <- get_coef_range_or()

if ( (TRUE == use_same_extra_term_coef_vec_all_sim) &
     (gen_extra_confounder_in_study | gen_extra_confounder_in_negctrl)) {
    print("> DEBUG: Generating coefficient of U1 outside of simulation loop.")

    beta_u_vec_template <- rep(NA_real_, times = num_outcomes)

    if (equi_confounding) {
        print(glue("> DEBUG: U1 is equi-confounding. Setting same effect size ",
                   "{signif(beta_u_for_all_outcomes, 3)} to all outcomes."))
        beta_u_vec_template <- rep(beta_u_for_all_outcomes, times = num_outcomes)

    } else {
        beta_u_vec_template <- set_coef_vals(
            initial_coef_vec = rep(NA_real_, times = num_outcomes),
            prop_negctrl_to_set = prop_beta_u_with_effect,
            val_range = beta_u_range)
    }
    if (debug_intrpt) { print("> DEBUG `beta_u_vec_template`") ; print(beta_u_vec_template) ; browser() }

} # End of code that generate coefficients of extra term *outside* of simulation loop.

####
### DataGen, Outcome model: Add to outcome model coef template: Quad term --------------------------
####

use_same_quad_term_coef_all_sim <- ifelse(
    !is.na(settings_file),
    yes = get_config("use_same_quad_term_coef_all_sim", default_val = FALSE) %>% as.logical(),
    no = FALSE)

static_quad_term_coef <- log(2.5)
if (outcome_model_coef_confounder_same_eff_all_outcomes_all_sims) {
    static_quad_term_coef <- confounder_to_all_outcomes_all_sims_eff
}

beta_x1sqr_vec_template <- rep(NA_real_, times = num_outcomes)
beta_x1sqr_range <- get_coef_range_or()

if (isTRUE(quad_term_sim_scenario)) {
    print("> DEBUG: Adding coefficients of quad term to OUTCOME model")
    beta_x1sqr_vec_template[seq(1L, num_outcomes)] <- 0.0

    beta_x1sqr_vec_template <- beta_x1sqr_vec_template %>%
        if_cond_do(use_same_quad_term_coef_all_sim, fun =  { function(coefvec) {
            rep(static_quad_term_coef, times = length(coefvec))
        }}) %>%
        if_cond_do(!use_same_quad_term_coef_all_sim & gen_study_with_x1sqr, fun = {
            function(coefvec) {
                # Outcome of interest depends on quad term.
                coefvec[[get_outcome_of_interest_study_num()]] <- log(
                runif(1, min = beta_x1sqr_range[['min_or']],
                         max = beta_x1sqr_range[['max_or']]))
                return(coefvec)  # End Y(X1Sqr)
            }
        }) %>%
        if_cond_do(!use_same_quad_term_coef_all_sim & gen_negctrls_with_x1sqr, fun = {
            function(coefvec) {
                # Negative control outcomes depend on quad term.
                coefvec[-get_outcome_of_interest_study_num()] <- log(
                    runif(num_negctrls, min = beta_x1sqr_range[['min_or']],
                          max = beta_x1sqr_range[['max_or']]))
                return(coefvec) # End W(X1Sqr)
            }
        })

    if (debug_intrpt) {
        print(glue("> DEBUG: Co-efficient template for quad term in OUTCOME model `beta_x1sqr_vec_template`
                        `use_same_quad_term_coef_all_sim`: {use_same_quad_term_coef_all_sim} ,
                        `gen_study_with_x1sqr`: {gen_study_with_x1sqr} ,
                        `gen_negctrls_with_x1sqr` : {gen_negctrls_with_x1sqr}"))
        print(beta_x1sqr_vec_template)
        browser()
    }
} # End of specifying outcome model reg coefficients for quad term.

####
### DataGen, Outcome model: Add to outcome model coef template: X1*X2 ------------------------------
####

use_same_x1x2_term_coef_all_sim <- ifelse(
    !is.na(settings_file),
    yes = get_config("use_same_x1x2_term_coef_all_sim", default_val = "FALSE") %>% as.logical(),
    no = FALSE)

static_x1x2_term_coef <- log(2.5)
if (outcome_model_coef_confounder_same_eff_all_outcomes_all_sims) {
    static_x1x2_term_coef <- confounder_to_all_outcomes_all_sims_eff
}

beta_x1x2_vec_template <- rep(NA_real_, times = num_outcomes)
beta_x1x2_range <- get_coef_range_or()

if (isTRUE(intr_twoconfounders_sim_scenario)) { # settings relevant to this scenario
    print("> DEBUG: Adding coefficients of X1*X2 to OUTCOME model TEMPLATE")

    beta_x1x2_vec_template[seq(1L, num_outcomes)] <- 0.0

    beta_x1x2_vec_template <- beta_x1x2_vec_template %>%
        # The same X1*X2 coef. in ALL simulations
        if_cond_do(use_same_x1x2_term_coef_all_sim, fun =  { function(coefvec) {
            rep(static_x1x2_term_coef, times = length(coefvec))
        }}) %>%
        # Generate random X1*X2 coef. in per-sim template
        if_cond_do( (!use_same_x1x2_term_coef_all_sim & gen_study_with_x1x2), fun = { function(coefvec) {
                # Outcome of interest
                coefvec[[get_outcome_of_interest_study_num()]] <- log(
                    runif(1, min = beta_x1x2_range[['min_or']],
                             max = beta_x1x2_range[['max_or']]))
                return(coefvec)
        }}) %>%
        if_cond_do(!use_same_x1x2_term_coef_all_sim & gen_negctrls_with_x1x2, fun = { function(coefvec) {
                # Negative controls
                coefvec[-get_outcome_of_interest_study_num()] <- log(
                    runif(num_negctrls,
                          min = beta_x1x2_range[['min_or']], max = beta_x1x2_range[['max_or']]) )
                return(coefvec)
        }})

    if (debug_intrpt) {
        print(glue("> DEBUG: Co-efficient template for X1*X2 term in OUTCOME model `beta_x1x2_vec_template`.
                            `use_same_x1x2_term_coef_all_sim` : {use_same_x1x2_term_coef_all_sim},
                            `gen_study_with_x1x2`: {gen_study_with_x1x2} ,
                            `gen_negctrls_with_x1x2` : {gen_negctrls_with_x1x2}"))
        print(beta_x1x2_vec_template)
        browser()
    }
} # End TRUE == intr_twoconfounders_sim_scenario

####
### DataGen, Outcome model template: Lack of positivity case ---------------------------------------
####

use_same_trtpos_deltavar_all_sim <- ifelse(
    !is.na(settings_file),
    yes = get_config("use_same_trtpos_deltavar_all_sim", default_val = FALSE) %>% as.logical(),
    no = FALSE)

beta_trtpos_deltavar_vec_template <- vector(mode = "numeric", length = num_outcomes)
if (lackpositivity_sim_scenario) {
    print("> DEBUG: Outcome model: Inflating coef of covar that crates lack of positivity")

    beta_trtpos_deltavar_vec_template <- rep(0.0, times = num_outcomes)
    beta_trtpos_deltavar_range <- get_coef_range_or()

    static_trtpos_deltavar_coef <- log(
        (max = beta_trtpos_deltavar_range[['max_or']]) * lackpositivity_coef_inflation_fct)

    beta_trtpos_deltavar_vec_template <- beta_trtpos_deltavar_vec_template %>%
        # The same coef. in ALL simulations
        if_cond_do(use_same_trtpos_deltavar_all_sim, fun =  { function(coefvec) {
            rep(static_trtpos_deltavar_coef, times = length(coefvec))
        }}) %>%
        if_cond_do(!use_same_trtpos_deltavar_all_sim & !gen_negctrls_with_trtpos,
            fun = { function(coefvec) {
                # Outcome of interest ONLY
                coefvec[[get_outcome_of_interest_study_num()]] <- log(
                    runif(1, min = beta_trtpos_deltavar_range[['max_or']] - 0.01,
                             max = beta_trtpos_deltavar_range[['max_or']]) * lackpositivity_coef_inflation_fct)
                return(coefvec)
            }
        }) %>%
        # Generate random inflated coef. in per-sim template
        if_cond_do(!use_same_trtpos_deltavar_all_sim & gen_negctrls_with_trtpos, fun = { function(coefvec) {
            # Outcome of interest and negative controls
            coefvec <- log(
                runif(length(coefvec),
                      min = beta_trtpos_deltavar_range[['max_or']] - 0.01,
                      max = beta_trtpos_deltavar_range[['max_or']]) *
                lackpositivity_coef_inflation_fct)
            return(coefvec)
            }
        })

    if (debug_intrpt) {
        print(glue("> DEBUG: Co-efficient template for covar affecting treatment assignment in
                      OUTCOME model `beta_trtpos_deltavar_vec_template`
                      `use_same_trtpos_deltavar_all_sim` : {use_same_trtpos_deltavar_all_sim},
                      `gen_negctrls_with_trtpos` : {gen_negctrls_with_trtpos}"))
        print(beta_trtpos_deltavar_vec_template)
        browser()
    }
} # End if(lackpositivity_sim_scenario)

####
### DataGen, Outcome model: Add to outcome model coef template: Measurement error -----
####

use_same_merror_term_coef_vec_all_sim <- ifelse(
    !is.na(settings_file),
    yes = get_config("use_same_merror_term_coef_vec_all_sim", default_val = FALSE) %>% as.logical(),
    no = FALSE)

beta_merror_var_vec_template <- vector(mode = "numeric", length = num_outcomes)
if (measurement_err_confounder_scenario) {
    print("> DEBUG: Outcome model: Inflating coef of measurement error covar")

    beta_merror_var_vec_template <- rep(0.0, times = num_outcomes)
    beta_merror_var_coef_range <- get_coef_range_or()

    static_merror_coef <- log(
        (max = beta_merror_var_coef_range[['max_or']]) * merror_coef_oddsratio_inflation_fct)

    beta_merror_var_vec_template <- beta_merror_var_vec_template %>%
        # The same coef. in ALL simulations
        if_cond_do(use_same_merror_term_coef_vec_all_sim, fun =  { function(coefvec) {
            rep(static_merror_coef, times = length(coefvec))
        }}) %>%
        if_cond_do(!use_same_merror_term_coef_vec_all_sim & !gen_negctrls_with_merror_covar,
            fun = { function(coefvec) {
                # Outcome of interest ONLY
                coefvec[[get_outcome_of_interest_study_num()]] <- log(
                    runif(1, min = beta_merror_var_coef_range[['max_or']] - 0.01,
                             max = beta_merror_var_coef_range[['max_or']]) * merror_coef_oddsratio_inflation_fct)
                return(coefvec)
            }
        }) %>%
        if_cond_do(
            !use_same_merror_term_coef_vec_all_sim & gen_negctrls_with_merror_covar,
            fun = { function(coefvec) {
                # Generate random inflated coef. in per-sim template
                # Outcome of interest and negative controls
                coefvec <- log(
                    runif(length(coefvec),
                          min = beta_merror_var_coef_range[['max_or']] - 0.01,
                          max = beta_merror_var_coef_range[['max_or']]) *
                                    merror_coef_oddsratio_inflation_fct)
                return(coefvec)
            }
        })
} # End adding measurement error covar `measurement_err_confounder_scenario`

####
### DataGen, Outcome model: Add to outcome model coef template: X1*Z -------------------------------
####

use_same_x1z_term_coef_all_sim <- FALSE
static_x1z_term_coef <- log(2.5)
if (outcome_model_coef_confounder_same_eff_all_outcomes_all_sims) {
    static_x1z_term_coef <- confounder_to_all_outcomes_all_sims_eff
}

if (isTRUE(intr_conftrt_sim_scenario)) { # settings relevant to this scenario
    print("> DEBUG: Adding coefficients of X1*Z to OUTCOME model")

    beta_x1z_vec_template <- rep(0.0, times = num_outcomes)
    beta_x1z_range <- get_coef_range_or()

    beta_x1z_vec_template <- beta_x1z_vec_template %>%
        # The same X1*Z coef. in ALL simulations
        if_cond_do(use_same_x1z_term_coef_all_sim, fun =  { function(coefvec) {
            rep(static_x1z_term_coef, times = length(coefvec))
        }}) %>%
        # Generate random X1*Z coef. in per-sim template
        if_cond_do(!use_same_x1z_term_coef_all_sim, fun = { function(coefvec) {
            # Outcome of interest
            coefvec[[get_outcome_of_interest_study_num()]] <- log(
                    runif(1, min = beta_x1z_range[['min_or']],
                             max = beta_x1z_range[['max_or']]))

            # Negative controls
            coefvec[-get_outcome_of_interest_study_num()] <- log(
                runif(num_negctrls,
                      min = beta_x1z_range[['min_or']], max = beta_x1z_range[['max_or']]) )
            return(coefvec)
            }
        })

    if (debug_intrpt) {
        print("> DEBUG: Co-efficient template for X1*Z term in OUTCOME model, `beta_x1z_vec_template`")
        print(beta_x1z_vec_template)
        browser()
    }
} # End TRUE == intr_conftrt_sim_scenario

##
# Simulation starts ------------------------------------------------------------
##

pkgs_dep <- c("dplyr", "tibble", "tidyselect", "vctrs", "purrr", "glue", "WeightIt", "EmpiricalCalibration")
Sys.setenv(MKL_NUM_THREADS = 1) # MKL
Sys.setenv(OMP_NUM_THREADS = 1) # OpenMP and MKL

if (n_sim > parallelise_sim_threshold) {
    cluster <- makeCluster(nprocesses)
    registerDoParallel(cluster)
} else {
   registerDoSEQ()
}

print("> Running simulation")
sim_start_time <- lubridate::now()
print(sim_start_time)

sim_out <- foreach(i = 1:n_sim, .packages = pkgs_dep,
    .verbose = TRUE) %dopar% {

    sim_number <- i

    # List to return from current simulation iteration.
    sim_output <- list()

    ###
    ## DataGen PER sim: Generates values of confounders --------------------------------------------
    ###

    print(glue("> DEBUG: Generating confounder values for with {num_negctrls} negative controls."))

    confounders_val_tbl <- gen_confounders_from_spec(n_obs = n_obs,
                                                 confounder_spec = confounder_spec)

    ### DataGen PER sim: Generates values for quad term -----
    ### This step is necessary because quad term is not part of `confounder_spec`
    ### because it's not an independent term.

    if (isTRUE(quad_term_sim_scenario)) {
        # Generate values for quadratic term here, since it depends on
        # an existing covariate (measured confounder).
        print(glue("> DEBUG: Sim #{sim_number}: Calculating values of quad term."))

        quadterm_basis_confounder_varname <- paste0(get_confounder_colname_prefix(), 1)
        quadterm_confounder_varname <- get_quadterm_colname()

        print(glue("> DEBUG: Generating data for quad term `{quadterm_confounder_varname}` ",
                   "using `{quadterm_basis_confounder_varname}`"))

        x1sqr_vals <- (confounders_val_tbl %>% pull(quadterm_basis_confounder_varname))^2
        confounders_val_tbl <- confounders_val_tbl %>% mutate(
                                 !!sym(quadterm_confounder_varname) := x1sqr_vals)

    } # End generating quad term data.

    ### DataGen PER sim: Generates values for X1*X2 interaction term -----
    ### Similar to the quadratic term case, this step is necessary because X1*X2 is not part
    ### of `confounder_spec`.
    if (isTRUE(intr_twoconfounders_sim_scenario)) {
        first_conf_name <- paste0(get_confounder_colname_prefix(), 1)
        second_conf_name <- paste0(get_confounder_colname_prefix(), 2)

        print(glue("> DEBUG: Sim #{sim_number}: Generating data for X1*X2 interaction term ",
                   "using {first_conf_name} and {second_conf_name}"))

        x1x2_intr_vals <- (confounders_val_tbl %>% pull(first_conf_name)) *
                          (confounders_val_tbl %>% pull(second_conf_name))

        confounders_val_tbl <- confounders_val_tbl %>% mutate(
            !!sym(get_x1x2intr_colname()) := x1x2_intr_vals)

    } # End generating VALUES for X1*X2 interaction, if necessary

    ### DataGen PER sim: Add measurement error -----
    ### This step is performed post-data generation using `confounders_spec`
    if (measurement_err_confounder_scenario & rand_merror_gaussian_params) {
        # Randomly creates mean and sd of measurement error
        print(glue("> DEBUG: Sim #{sim_number}: Measurement error, selecting params for rand error"))
        merror_gaussian_mean <- runif(1, min = rand_merror_gaussian_mean_min,
                                         max = rand_merror_gaussian_mean_max)
        rand_merror_gaussian_mean_deviation_pct <- runif(1, min = rand_merror_gaussian_mean_deviation_min_pct,
                                            max = rand_merror_gaussian_mean_deviation_max_pct)
        rand_merror_gaussian_mean_deviation <- merror_gaussian_mean *
                                                (rand_merror_gaussian_mean_deviation_pct / 100)
        merror_gaussian_sd <- calc_gaussian_sd_from_mean_deviation(merror_gaussian_mean,
                                                rand_merror_gaussian_mean_deviation)

        merror_params <- tibble(
            gaussian_mean = merror_gaussian_mean,
            deviation_from_mean_pct = rand_merror_gaussian_mean_deviation_pct,
            gaussian_sd = merror_gaussian_sd)
        sim_output[["merror_gaussian_params"]] <- merror_params

        print(glue("> DEBUG: Sim #{sim_number}: Measurement error. Gaussian error
                                                mean: {signif(merror_gaussian_mean, 4)},
                                                  sd: {signif(merror_gaussian_sd, 4)}"))
    }

    if (measurement_err_confounder_scenario) {
        print(glue("> DEBUG: Sim #{sim_number}: Adding measurement error"))

        covar_vec <- confounders_val_tbl %>% pull(get_merror_confounder_colname())
        merr_cover_vec <- covar_vec + rnorm(length(covar_vec),
            mean = merror_gaussian_mean, sd = merror_gaussian_sd)

        confounders_val_tbl <- confounders_val_tbl %>% mutate(
            !!sym(get_merror_confounder_colname()) := merr_cover_vec)
    } # End adding measurement error of confounder

    if (debug_intrpt) { browser() }

    ###
    ## DataGen PER sim: Generate values in TREATMENT MODEL ----------------------------------------
    ###

    ### DataGen PER sim: Creates the coefficients for treatment model per simulation ----
    trt_model_coefs_range <- get_coef_range_or()

    trt_model_coefs_vec <- vector("numeric", length = ncol(confounders_val_tbl))

    if (all_sim_use_trt_model_template) {
        # Uses the template across all simulations
        trt_model_coefs_vec <- trt_model_coefs_template_vec

    } else {
        print(glue("> DEBUG: Sim #{sim_number}: DataGen Generating simulation specific treatment reg coefs"))
        trt_model_coefs_vec <- confounder_spec %>% gen_trt_model_regression_coefs(
            coef_range = trt_model_coefs_range)

        # Build the treatment model coef vector: only for those terms NOT in confounder spec.
        trt_model_coefs_vec <- trt_model_coefs_vec %>%
            # Adding coefficient for quadratic term (since quad term is not in confounder spec
            if_cond_do(condition = isTRUE(quad_term_sim_scenario) & isTRUE(x1sqr_in_trt_model),
                fun = { function(coef_vec) {
                     print(glue("> DEBUG: Sim #{sim_number}: Generating ",
                        "reg. co-efficient for quadratic term in TREATMENT model."))

                     quadterm_coef_trtmodel <- log(
                         runif(1, min = trt_model_coefs_range[['min_or']],
                                  max = trt_model_coefs_range[['max_or']]))

                     coef_vec <- vctrs::vec_c(coef_vec,
                                              !!sym(get_quadterm_colname()) := c(quadterm_coef_trtmodel))
                     return(coef_vec)
                    }
            }) %>%
            if_cond_do(condition = isTRUE(intr_twoconfounders_sim_scenario), fun = { function(coef_vec) {
                print(glue("> DEBUG: Sim #{sim_number}: Generating ",
                                     "per-sim reg. co-efficient for X1*X2 in TREATMENT model."))

                coef_val <- log(runif(1, min = trt_model_coefs_range[['min_or']],
                                         max = trt_model_coefs_range[['max_or']]))

                coef_vec <- vctrs::vec_c(coef_vec, !!sym(get_x1x2intr_colname()) := c(coef_val))
                return(coef_vec)
            }}) %>%
            if_cond_do(condition = lackpositivity_sim_scenario, fun = { function(coef_vec) {
                print(glue("> DEBUG: Sim #{sim_number}: Generating ",
                           "coef for covar that induces lack of positivity in TREATMENT model."))

                coef_val <- log(trt_model_coefs_range[['max_or']] * lackpositivity_coef_inflation_fct)
                # Replacing existing value because this term is part of confounder spec
                coef_vec[[get_trtpos_deltavar_colname()]] <- c(coef_val)
                return(coef_vec)
            }})

    } # End ifelse(all_sim_use_trt_model_template)

    # Returns reg coefs for outcome model if required
    if (persist_sim_out_data) {
        sim_output[['trt_model_coef']] <- as.list(trt_model_coefs_vec) %>% as_tibble()
    }

    print(glue("> DEBUG: Sim #{sim_number}: Generating values for treatment"))

    ### Sets up confounders data matrix
    confounders_val_tbl_trtmodel <- confounders_val_tbl

    confounders_val_tbl_trtmodel <- confounders_val_tbl_trtmodel %>%
        # Conditional treatment model for quad term scenario — un-necessary
        if_cond_do(condition = ((TRUE == quad_term_sim_scenario) & (FALSE == x1sqr_in_trt_model)),
                   fun = { function(design_mat) {
                print(glue("> DEBUG: Sim #{sim_number}: ",
                    "`x1sqr_in_trt_model` is {x1sqr_in_trt_model}. ",
                    "In TREATMENT model, remove `X1Sqr` from design matrix."))
                confounders_val_tbl_trtmodel %>% select(-.data[[get_quadterm_colname()]])
            }
        })

    ##### Adds intercept term to confounder data
    confounder_design_mat <- confounders_val_tbl_trtmodel %>%
        mutate(INTR = rep(1.0, times = n_obs)) %>% relocate(INTR) %>%
        as.matrix()

    ### Checking design matrix and coef vector of treatment model.

    print(glue("> DEBUG: Sim #{sim_number}: Checking design matrix and vector of coef compatibility"))
    #print(head(confounder_design_mat, 5))
    #print(trt_model_coefs_vec)

    trt_designmat_compat_coef_vec <- assertthat::are_equal(ncol(confounder_design_mat), length(trt_model_coefs_vec))
    if (!isTRUE(trt_designmat_compat_coef_vec)) {
        ncol_design_mat <- ncol(confounder_design_mat)
        length_coef_vec <- length(trt_model_coefs_vec)
        err_msg <- glue("!! ERROR: Sim #{sim_number}: TREATMENT model design matrix and vector of coef ",
                        "have different lengths; not `comfortable`. ",
                        "Design matrix has {ncol_design_mat} columns, length of coef vector is {length_coef_vec}.")
        stop(err_msg)
    }

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: Design matrix for TREATMENT model"))
        print(head(confounder_design_mat))
        print(glue("> DEBUG: Sim #{sim_number}: Vector of regression coefs for TREATMENT model"))
        print(trt_model_coefs_vec)
        #print(glue("> DEBUG: Sim #{sim_number}: TREATMENT model design matrix `comfortable` ",
        #           "with coef vector: {trt_designmat_compat_coef_vec}"))
        browser()
    }

    # Prob. of treated as a combination of covariates.
    trt_model_eta <- confounder_design_mat %*% trt_model_coefs_vec
    pr_trt <- 1 / (1 + exp(-1 * trt_model_eta))
    Z <- rbinom(n_obs, size = 1L, prob = pr_trt)

    if (lackpositivity_sim_scenario) {
        #browser()
        # Manually assign treatment based on value of covariate.
        print(glue("> DEBUG: Sim #{sim_number}: Generating treatment status based ",
                   "on covar value to simulate lack of positivity"))

        # Using setup similar to Kang J, Chan W, Kim MO, Steiner P, 2016.
        # DOI 10.5351/CSAM.2016.23.1.001
        # Directly manipulates the propensity score vector.
        x_trt_vec <- confounder_design_mat[, get_trtpos_deltavar_colname()]

        x_trt_0_indices  <- x_trt_vec <= -0.8
        x_trt_1_indices  <- x_trt_vec >= 0.5
        x_trt_w_indicies <- (x_trt_vec > -0.8) & (x_trt_vec < 0.5)

        pscore_vec <- rep(NA_real_, n_obs)
        pscore_vec[x_trt_0_indices] <- 0.0
        pscore_vec[x_trt_1_indices] <- 1.0

        x_trt_w_eta <- confounder_design_mat[x_trt_w_indicies, ] %*% trt_model_coefs_vec
        pscore_vec[x_trt_w_indicies] <- x_trt_w_eta

        Z <- 1 * (pscore_vec > runif(n_obs))

    } # End generating treatment for non-positivity

    ### DataGen PER sim: After gen. of trt values, generates values for X1*Z term -----
    if (isTRUE(intr_conftrt_sim_scenario)) {
        print(glue("> DEBUG: Sim #{sim_number}: Calculating values of X1*Z term."))

        intr_conftrt_basis_confounder_varname <- paste0(get_confounder_colname_prefix(), 1)
        trt_varname <- get_trt_colname()

        print(glue("> DEBUG: Generating data for X1*Z term`"))
        x1z_vals <- (confounders_val_tbl %>% pull(intr_conftrt_basis_confounder_varname)) * Z

        confounders_val_tbl <- confounders_val_tbl %>%
            mutate(!!sym(get_x1_zintr_colname()) := x1z_vals)
    } # End TRUE == intr_conftrt_sim_scenario

    # Adds treatment values to confounder table
    trt_data <- confounders_val_tbl %>% mutate(Z = Z)

    # Change the column order so that it's in the same order for outcome model.
    trt_data <- trt_data %>% relocate(Z,
        .before = paste0(get_confounder_colname_prefix(), 1L))

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: `trt_data`"))
        print(head(trt_data, 5)) ; browser()
    }

    pscore_terms <- colnames(confounders_val_tbl) %>%
        if_cond_do(condition = ((TRUE == quad_term_sim_scenario) & (FALSE == x1sqr_in_trt_model)),
                   fun = { function(covarnames) {
                       colname_is_quadterm_vec <- covarnames %in% get_quadterm_colname()
                       covarnames[!colname_is_quadterm_vec]
                   }
        })
    pscore_formula_rhs <- paste0(pscore_terms, collapse = "+")
    pscore_formula <- paste("Z", pscore_formula_rhs, sep = " ~ ")
    print(glue("> DEBUG: Sim #{sim_number}: Treatment model Calculating IPW weights {pscore_formula}"))
    stablised_pscores_untrimmed_obj <- weightit(as.formula(pscore_formula), data = trt_data, method = "ps", stabilize = TRUE)

    # Trim weights at 99 quantile
    stablised_pscores_obj <- WeightIt::trim(stablised_pscores_untrimmed_obj, at = 0.99)

    if (debug_intrpt) {
        # Balance diagnostic graphs
        trt_pscore_plot <- cobalt::bal.plot(stablised_pscores_obj, which = "both")
        print(trt_pscore_plot) ; browser()
    }

    #sim_output[["trt_pscores"]] <- list(tibble(sim_num = rep(sim_number, times = length(Z)), Z = Z, pscores = pscores_obj$ps))
    sim_output[["trt_stablised_pscores"]] <- list(
        tibble(sim_num = rep(sim_number, times = length(Z)),
            Z = Z, sweights = stablised_pscores_obj$weights))
    sim_output[['trt_assignment']] <- tibble(
        sim_num = sim_number, num_trt = sum(Z), num_not_trt = n_obs - sum(Z))

    ###
    ### DataGen PER sim Outcome model: Populate reg coef template with values ------------------
    ###

    outcome_model_coef_tbl <- outcome_model_coef_tbl_template
    true_trt_eff_persim <- trt_eff

    # Init parameters
    #same_all_measured_eff_to_all_outcomes_all_sim <- NA
    #same_measured_eff_to_all_outcomes_per_sim_iter <- NA

    # Each simulation iteration generates its own outcome model
    if (FALSE == use_outcome_model_coef_template_all_sims) {

        # In each simulation,
        # coef. of ALL measured confounders, for ALL outcomes, are the same.
        #same_all_measured_eff_to_all_outcomes_all_sim <- get_config(
        # "same_all_measured_eff_to_all_outcomes_all_sim", default_val = FALSE) %>% as.logical()
        sim_output[['same_all_measured_eff_to_all_outcomes_all_sim']] <- same_all_measured_eff_to_all_outcomes_all_sim

        common_coef_all_confounders_all_outcomes <- log(
            runif(1L, min = measured_conf_coef_range[['min_or']],
                max = measured_conf_coef_range[['max_or']]))
        sim_output[['common_coef_all_confounders_all_outcomes']] <- common_coef_all_confounders_all_outcomes

        # In each simulation,
        # each measured confounders have the same effect size across all outcomes,
        # BUT the coef. differ across confounders
        # (unlike `same_all_measured_eff_to_all_outcomes_all_sim`)
        #same_measured_eff_to_all_outcomes_per_sim_iter <- get_config(
        #   "same_measured_eff_to_all_outcomes_per_sim_iter", default_val = FALSE) %>% as.logical()
        sim_output[['same_measured_eff_to_all_outcomes_per_sim_iter']] <- same_measured_eff_to_all_outcomes_per_sim_iter

        print(glue("> DEBUG: Sim #{sim_number}: Generating regression coefficients for MEASURED confounder X in OUTCOME MODEL."))
        print(glue("> DEBUG: Sim #{sim_number}: MEASURED confounders coef same for ALL outcomes, across ALL sims: {same_all_measured_eff_to_all_outcomes_all_sim}"))
        print(glue("> DEBUG: Sim #{sim_number}: MEASURED confounders coef same for ALL outcomes, PER sim: {same_measured_eff_to_all_outcomes_per_sim_iter}"))

        # Generates per-simulation 'true' treatment effect
        true_trt_eff_persim <- log( runif(1, min = trt_eff_range[['min_or']],
                                             max = trt_eff_range[['max_or']]))
        trt_eff_persim_vec <- c(true_trt_eff_persim,
            rep(0.0, times = nrow(outcome_model_coef_tbl) - 1L) )

        outcome_model_trt_colname <- paste0(outcome_model_coef_prefix,
            get_trt_colname())
        outcome_model_coef_tbl <- outcome_model_coef_tbl %>%
                mutate(!!outcome_model_trt_colname := trt_eff_persim_vec)

        outcome_model_measured_confounders_names <- paste0(outcome_model_coef_prefix,
            seq.int(0, to = num_notpriout_mconf + 1L))

        #browser()

        for (col_index in seq_along(outcome_model_measured_confounders_names)) {
            # Loop over each measured confounder.

            measured_conf_colname <-  outcome_model_measured_confounders_names[[col_index]]

            current_measured_conf_beta_vec <- outcome_model_coef_tbl %>%
                pull(.data[[measured_conf_colname]])

            beta_x_vec <- rep(NA_real_, times = length(current_measured_conf_beta_vec))

            if (same_all_measured_eff_to_all_outcomes_all_sim) {
                 print(glue("> DEBUG: Sim #{sim_number}: DataGen ",
                     "All confounder-outcome coef are the same."))

                beta_x_vec <- set_coef_vals(initial_coef_vec = current_measured_conf_beta_vec,
                                            prop_negctrl_to_set =
                                            1.0, val = common_coef_all_confounders_all_outcomes)

            } else if (same_measured_eff_to_all_outcomes_per_sim_iter) {

                print(glue("> DEBUG: Sim #{sim_number}: DataGen ",
                    "Each measured confounder-outcome pair has same coef., ",
                    "but DIFFER across confounders."))
                #print(glue("> DEBUG: Sim #{sim_number}: DataGen ",
                #    "Generating coef. of {measured_conf_colname} to ALL outcomes"))

                coef_to_outcomes <- log(
                    runif(1L, min = measured_conf_coef_range[['min_or']],
                              max = measured_conf_coef_range[['max_or']]))

                beta_x_vec <- set_coef_vals(initial_coef_vec = current_measured_conf_beta_vec,
                                            prop_negctrl_to_set = 1.0,
                                            val = coef_to_outcomes)

                if (length(unique(beta_x_vec)) != 1) {
                    stop(glue("!! ERROR: Sim #{sim_number}: Coefs for {measured_conf_colname} not all the same, ",
                              "`same_measured_eff_to_all_outcomes_per_sim_iter` : {same_measured_eff_to_all_outcomes_per_sim_iter}"))
                }

            } else {

                beta_x_vec <- set_coef_vals(initial_coef_vec = current_measured_conf_beta_vec,
                                            prop_negctrl_to_set = 1.0,
                                            val_range = measured_conf_coef_range)
            }

            #print(glue("> DEBUG: Sim #{sim_number}: Coef for {measured_conf_colname}"))
            #print(beta_x_vec)

            outcome_model_coef_tbl <- outcome_model_coef_tbl %>%
                mutate(!!measured_conf_colname := beta_x_vec)

        } # Looping over non-treatment measured confounders
    } # End if (FALSE == use_outcome_model_coef_template_all_sims) , i.e., no template to use.

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: INITIAL outcome model reg coef.
                    `use_outcome_model_coef_template_all_sims` : {use_outcome_model_coef_template_all_sims}
                    `same_measured_eff_to_all_outcomes_per_sim_iter` : {same_measured_eff_to_all_outcomes_per_sim_iter}
                   "))
        print(outcome_model_coef_tbl) ; browser()
    }

    ###
    ### DataGen PER sim Outcome model: Reg coef for potential unmeasured confounder U1 ---------
    ###

    if (gen_extra_confounder_in_study | gen_extra_confounder_in_negctrl) {
        # Generate coefficients for potential unmeasured confounder

        print(glue("> DEBUG: Sim #{sim_number}: Setting up regression coefficients for potential unmeasured confounder U in OUTCOME MODEL."))
        print(glue("> DEBUG: Sim #{sim_number}: Proportional of negative controls affected by U: {signif(prop_beta_u_with_effect, 3)}, ",
                   "which is {(num_negctrls * prop_beta_u_with_effect) %>% ceiling() %>% as.integer()} neg ctrls."))

        # Adds initial column of empty values for potential unmeasured confounder, U
        outcome_model_coef_tbl <- outcome_model_coef_tbl %>%
            add_one_outcome_model_coef(coef_name = potential_unmeasured_confounder_varname)

        beta_u_vec <- vector("numeric", length = num_outcomes)
        if (FALSE == use_same_extra_term_coef_vec_all_sim) {
            print(glue("> DEBUG: Sim #{sim_number}: GENERATING regression coefficients for ",
                       "potential unmeasured confounder U in OUTCOME MODEL."))

            if (TRUE == equi_confounding) {
                print(glue("> DEBUG: Sim #{sim_number}: Equi-confounding. Same beta_u, ",
                           "effect of pontential unmeasured confounder, for all outcomes."))
                beta_u_vec <- set_coef_vals(initial_coef_vec = outcome_model_coef_tbl %>%
                                                pull(.data[[potential_unmeasured_confounder_varname]]),
                                            prop_negctrl_to_set = 1.0,
                                            val = beta_u_for_all_outcomes)

            } else {
                # Effect of potential unmeasured confounder U, initially as randomly set of values
                beta_u_vec <- set_coef_vals(
                    initial_coef_vec = outcome_model_coef_tbl %>% pull(.data[[potential_unmeasured_confounder_varname]]),
                    prop_negctrl_to_set = prop_beta_u_with_effect,
                    val_range = beta_u_range)
            }

        } else if (use_same_extra_term_coef_vec_all_sim &
            same_measured_eff_to_all_outcomes_per_sim_iter) {

            print(glue("> DEBUG: Sim #{sim_number}: DataGen ",
                "Coef for potential unmeasured confounder U ",
                "same for ALL outcome per simulation. Equi-confounding by default."))

            beta_u_for_all_outcomes_current_sim <- log(
                runif(1, min = beta_u_range[['min_or']], max = beta_u_range[['max_or']]) )

            beta_u_vec <- set_coef_vals(
                initial_coef_vec = outcome_model_coef_tbl %>%
                                   pull(.data[[potential_unmeasured_confounder_varname]]),
                prop_negctrl_to_set = 1.0,
                val = beta_u_for_all_outcomes_current_sim)

        } else {
            beta_u_vec <- beta_u_vec_template

        } # if FALSE == use_same_extra_term_coef_vec_all_sim

        ## Outcome model: Negative controls are not affected by unmeasured confounder U
        if (!gen_extra_confounder_in_negctrl) {
            print(glue("> DEBUG: Sim #{sim_number} ",
                       "`gen_extra_confounder_in_negctrl:` {gen_extra_confounder_in_negctrl} ",
                       "Removing association between potential unmeasured confounder U and negative controls."))

            beta_u_vec[2:length(beta_u_vec)] <- 0.0

        } # End if !gen_extra_confounder_in_negctrl

        outcome_model_coef_tbl <- outcome_model_coef_tbl %>%
            mutate(!!potential_unmeasured_confounder_varname := beta_u_vec)

        #if (debug_intrpt) { print(outcome_model_coef_tbl) ; browser() }
    } # End of specifying outcome model reg coefficients for unmeasured confounding case.

    ###
    ### DataGen PER sim Outcome model: Reg coef for quad term ----------------------------------
    ###
    if (isTRUE(quad_term_sim_scenario)) {
        print(glue("> DEBUG: Sim #{sim_number}: Setting up regression coefficients for ",
                "quadratic term in OUTCOME MODEL."))

        # By default, use whatever value that is in template
        beta_x1sqr_coef_vec <- beta_x1sqr_vec_template

        if (FALSE == use_same_quad_term_coef_all_sim) {
            # Randomly sample coefficients for all outcomes, in each sim iteration.

            if (isTRUE(gen_study_with_x1sqr)) {
                # Outcome of interest depends on quad term.
                beta_x1sqr_coef_vec[[get_outcome_of_interest_study_num()]] <-
                    log(runif(1, min = beta_x1sqr_range[['min_or']], max = beta_x1sqr_range[['max_or']]))
            } # End Y(X1Sqr)

            if (isTRUE(gen_negctrls_with_x1sqr)) {
                # Negative control outcomes depend on quad term.
                beta_x1sqr_coef_vec[-get_outcome_of_interest_study_num()] <-
                    log(runif(num_negctrls, min = beta_x1sqr_range[['min_or']], max = beta_x1sqr_range[['max_or']]))
            } # End W(X1Sqr)
        }
        else if (use_same_quad_term_coef_all_sim & same_measured_eff_to_all_outcomes_per_sim_iter) {
            # Same coef. to all outcome in each sim. iteration.
            # Ignores `gen_study_with_x1sqr` and `gen_negctrls_with_x1sqr`
            print(glue("> DEBUG: Sim #{sim_number}: DataGen ",
                       "Coef of quad term X1Sqr same for ALL outcome, in each simulation."))

            coef_for_all_outcomes <- log( runif(1,
                                                min = beta_x1sqr_range[['min_or']],
                                                max = beta_x1sqr_range[['max_or']]))

            beta_x1sqr_new <- rep(NA, times = num_outcomes) # Allows `set_coef_vals()` to set coef for outcome of interest
            beta_x1sqr_coef_vec <- set_coef_vals(
                initial_coef_vec = beta_x1sqr_new,
                prop_negctrl_to_set = 1.0,
                val = coef_for_all_outcomes)
        }

        outcome_model_coef_tbl <- outcome_model_coef_tbl %>%
            mutate(!!get_quadterm_colname() := beta_x1sqr_coef_vec)

    } # End of specifying outcome model reg coefficients for quad term.

    ###
    ### DataGen PER sim Outcome model: Reg coef for X1*X2 interaction -------------------------------
    ###
    if (isTRUE(intr_twoconfounders_sim_scenario)) {
        print(glue("> DEBUG: Sim #{sim_number}: Setting up per-sim regression coefficients for ",
                   "X1*X2 interaction term in OUTCOME MODEL."))

        # By default, use whatever value that is in template
        beta_x1x2_vec <- beta_x1x2_vec_template

        if (FALSE == use_same_x1x2_term_coef_all_sim) {
            # Generating a new set of coef. for the current simulation iteration.
            beta_x1x2_vec <- gen_outcome_model_coef_persim(
                coef_in_primaryoutcome = gen_study_with_x1x2,
                coef_in_negctrls = gen_negctrls_with_x1x2,
                coef_range = beta_x1x2_range,
                num_outcomes = num_outcomes,
                primary_outcome_num = 1L,
                coef_noeff_val = 0.0)

        } else {
            # Same coefficient for all outcomes
            beta_x1x2_vec <- gen_same_outcome_model_coef_persim(
                coef_range = beta_x1x2_range, num_outcomes = num_outcomes, coef_noeff_val = 0.0)
        }

        outcome_model_coef_tbl <- outcome_model_coef_tbl %>%
            mutate(!!get_x1x2intr_colname() := beta_x1x2_vec)

    } # End of specifying outcome model reg coefficients for X1*X2 interaction term.

    ###
    ### DataGen PER sim Outcome model: Reg coef for non-positivity var ---------
    ###
    if (isTRUE(lackpositivity_sim_scenario)) {
        print(glue("> DEBUG: Sim #{sim_number}: Setting up per-sim regression ",
                   "coefficients for non-positivity covar in OUTCOME MODEL."))

        # By default, use whatever value that is in template
        beta_xpostrt_vec <- beta_trtpos_deltavar_vec_template

        # Generates per-sim coef. is necessary
        if (FALSE == use_same_trtpos_deltavar_all_sim) {
            # Generating a new set of coef. for the current simulation iteration.
            coef_max_or <- beta_xpos_deltavar_range[['max_or']]
            new_coef_min_or <- coef_max_or - 0.1
            new_coef_max_or <- coef_max_or * lackpositivity_coef_inflation_fct

            if (gen_negctrls_with_trtpos) {
                beta_xpostrt_vec <- log(
                    runif(length(beta_trtpos_deltavar_vec_template),
                      min = new_coef_min_or, max = new_coef_max_or))
            } else {
                beta_xpostrt_vec[[get_outcome_of_interest_study_num()]] <- log(
                    runif(1L, min = new_coef_min_or, max = new_coef_max_or))
            }

        } else if (use_same_trtpos_deltavar_all_sim & same_measured_eff_to_all_outcomes_per_sim_iter) {
            # Same coef. to all outcome in each sim. iteration.
            beta_trtpos_deltavar_range_persim <- list()
            beta_trtpos_deltavar_range_persim[['min_or']] <- beta_trtpos_deltavar_range[['max_or']] - 0.01
            beta_trtpos_deltavar_range_persim[['max_or']] <- beta_trtpos_deltavar_range[['max_or']] * lackpositivity_coef_inflation_fct

            beta_xpostrt_vec <- gen_same_outcome_model_coef_persim(
                coef_range = beta_trtpos_deltavar_range_persim, num_outcomes = num_outcomes, coef_noeff_val = 0.0)
        }

        outcome_model_coef_tbl <- outcome_model_coef_tbl %>%
            mutate(!!get_trtpos_deltavar_colname() := beta_xpostrt_vec)

    } # End specify reg coef for lack of positivity var.

    ###
    ### DataGen PER sim Outcome model: Reg coef for measurement error term -----
    ###
    if (isTRUE(measurement_err_confounder_scenario)) {
         print(glue("> DEBUG: Sim #{sim_number}: Setting up per-sim regression ",
                   "coefficients for measurement error covar in OUTCOME MODEL."))

        # By default, use whatever value that is in template
        beta_merror_vec <- beta_merror_var_vec_template

        # Generates per-sim coef. is necessary
        if (FALSE == use_same_merror_term_coef_vec_all_sim) {
            # Generating a new set of coef. for the current simulation iteration.
            coef_max_or <- beta_merror_var_coef_range[['max_or']]
            new_coef_min_or <- coef_max_or - 0.1
            new_coef_max_or <- coef_max_or * merror_coef_oddsratio_inflation_fct

            if (gen_negctrls_with_merror_covar) {
                beta_merror_vec <- log(
                    runif(length(beta_merror_var_vec_template),
                      min = new_coef_min_or, max = new_coef_max_or))
            } else {
                beta_merror_vec[[get_outcome_of_interest_study_num()]] <- log(
                    runif(1L, min = new_coef_min_or, max = new_coef_max_or))
            }
        } else if (use_same_merror_term_coef_vec_all_sim & same_measured_eff_to_all_outcomes_per_sim_iter) {
            # Same coef. to all outcome in each sim. iteration.
            beta_range_persim <- list()
            beta_range_persim[['min_or']] <- beta_merror_var_coef_range[['max_or']] - 0.01
            beta_range_persim[['max_or']] <- beta_merror_var_coef_range[['max_or']] * merror_coef_oddsratio_inflation_fct

            beta_merror_vec <- gen_same_outcome_model_coef_persim(
                coef_range = beta_range_persim,
                num_outcomes = num_outcomes, coef_noeff_val = 0.0)
        }

        outcome_model_coef_tbl <- outcome_model_coef_tbl %>%
            mutate(!!get_merror_confounder_colname() := beta_merror_vec)

    } # End specifying per-sim reg. coef. for measurement error

    ###
    ### DataGen PER sim Outcome model: Reg coef for X1*Z interaction ---------
    ###
    if (isTRUE(intr_conftrt_sim_scenario)) {
        print(glue("> DEBUG: Sim #{sim_number}: Setting up regression coefficients for ",
                   "X1*Z term in OUTCOME MODEL."))
        # By default, use whatever value that is in template
        beta_x1z_coef_vec <- beta_x1z_vec_template

        if (FALSE == use_same_x1z_term_coef_all_sim) {
            # Randomly sample coefficients for all outcomes, in each sim iteration.

            # Outcome of interest
            beta_x1z_coef_vec[[get_outcome_of_interest_study_num()]] <-
                log(runif(1, min = beta_x1z_range[['min_or']], max = beta_x1z_range[['max_or']]))
            # Negative control outcomes depend on quad term.
            beta_x1z_coef_vec[-get_outcome_of_interest_study_num()] <-
                log(runif(num_negctrls, min = beta_x1z_range[['min_or']], max = beta_x1z_range[['max_or']]))

        }  else if (use_same_x1z_term_coef_all_sim & same_measured_eff_to_all_outcomes_per_sim_iter) {
            # Same coef. to all outcome in each sim. iteration.
            print(glue("> DEBUG: Sim #{sim_number}: DataGen ",
                       "Coef of X1*Z same for ALL outcome, in each simulation."))

            coef_for_all_outcomes <- log( runif(1,
                                                min = beta_x1z_range[['min_or']],
                                                max = beta_x1z_range[['max_or']]))

            beta_x1z_new <- rep(NA, times = num_outcomes) # Allows `set_coef_vals()` to set coef for outcome of interest
            beta_x1z_coef_vec <- set_coef_vals(
                initial_coef_vec = beta_x1z_new,
                prop_negctrl_to_set = 1.0,
                val = coef_for_all_outcomes)
        }

        outcome_model_coef_tbl <- outcome_model_coef_tbl %>%
            mutate(!!get_x1_zintr_colname() := beta_x1z_coef_vec)
    } # End of specifying per-sim outcome model reg coefficients for X1*Z term.

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: Final outcome model coef `outcome_model_coef_tbl`"))
        TRUE %>% if_cond_do(condition = gen_extra_confounder_in_study,
            fun = { function(cond) {
                print(glue("> DEBUG: Sim #{sim_number}: `gen_extra_confounder_in_study` {gen_extra_confounder_in_study}"))
                return(TRUE)
                }
            }) %>%
        if_cond_do(condition = gen_extra_confounder_in_negctrl,
            fun = { function(cond) {
                print(glue("> DEBUG: Sim #{sim_number}: `gen_extra_confounder_in_negctrl` {gen_extra_confounder_in_negctrl}"))
                return(TRUE)
                }
            }) %>%
        if_cond_do(condition = (quad_term_sim_scenario & gen_study_with_x1sqr),
            fun = { function(cond) {
                print(glue("> DEBUG: Sim #{sim_number}: `gen_study_with_x1sqr` {gen_study_with_x1sqr}"))
                return(TRUE)
                }
            }) %>%
        if_cond_do(condition = (quad_term_sim_scenario & gen_study_with_x1sqr),
            fun = { function(cond) {
                print(glue("> DEBUG: Sim #{sim_number}: `gen_negctrls_with_x1sqr` {gen_negctrls_with_x1sqr}"))
                return(TRUE)
                }
            })

        print(outcome_model_coef_tbl)
        browser()
    }

    # Returns reg coefs for outcome model if required
    if (persist_sim_out_data) {
        sim_output[['outcome_model_coef']] <- outcome_model_coef_tbl
    }

    ##
    # Outcome model: Generate outcomes ---------------------------------------------
    ##

    print(glue("> DEBUG: Sim #{sim_number}: Generating outcome values"))

    trt_data_design_matrix <- trt_data %>%
        mutate(INTR = rep(1L, times = nrow(trt_data))) %>%
        relocate("INTR") %>% as.matrix()

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: outcome model design matrix `trt_data_design_matrix`:"))
        print(head(trt_data_design_matrix, 5))
        print("")
        print(glue("> DEBUG: Sim #{sim_number}: outcome model coef `outcome_model_coef_tbl`:"))
        print(outcome_model_coef_tbl)
        browser()
    }

    outcomes_mat <- outcome_model_coef_tbl %>% apply(MARGIN = 1, FUN = function(row) {
        outcome_num <- row[[1]]

        # Removes the outcome ID
        coef_vec <- row[-1]

        outcome_model_eta <- trt_data_design_matrix %*% coef_vec
        pr_y1 <- 1 / (1 + exp(-1 * outcome_model_eta))
        Y <- rbinom(n_obs, size = 1L, prob = pr_y1)
    })
    colnames(outcomes_mat) <- paste0(get_outcome_colname(), seq.int(ncol(outcomes_mat)))

    if (debug_intrpt) { print(head(outcomes_mat)) ; browser() }

    ## Verify outcome model for outcome of interest
    print(glue("> DEBUG: Sim #{sim_number}: Verifying outcome of interest"))
    outcome_of_interest_data <- cbind(
        outcomes_mat[, 1], trt_data) %>% as_tibble()
    colnames(outcome_of_interest_data) <- c(
        get_outcome_colname(), colnames(trt_data))

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: Outcome of interest data `outcome_of_interest_data`"))
        print(head(outcome_of_interest_data, 5))
        browser()
    }

    if (persist_sim_out_data) {
        sim_output[['primary_outcome_data']] <- outcome_of_interest_data
    }

    trt_coef_verified <- verify_binary_outcome(data = outcome_of_interest_data,
        trt_var = get_trt_colname(), out_var = get_outcome_colname(),
        expected_trt_eff = trt_eff)

    print(glue("> DEBUG: Sim #{sim_number}: Treatment effect in outcome of interest matches {trt_coef_verified}"))

    ## Verify outcome model for a negative control
    if (num_negctrls > 0) {
        print(glue("> DEBUG: Sim #{sim_number}: Verifying first negative control"))

        negctrl_data_preproc <- cbind(
            outcomes_mat[, 2], trt_data) %>% as_tibble()
        colnames(negctrl_data_preproc) <- c(
            get_outcome_colname(), colnames(trt_data))

        trt_coef_verified <- verify_binary_outcome(data = negctrl_data_preproc,
            trt_var = get_trt_colname(), out_var = get_outcome_colname(),
            expected_trt_eff = 0.0)

        print(glue("> DEBUG: Sim #{sim_number}: Treatment effect in first negative outcome {trt_coef_verified}"))

        ### Also check for zero association with potential unmeasured confounder U1
        if (TRUE == gen_extra_confounder_in_negctrl) {
            negctrl_outcome_u_verified <- verify_binary_outcome(data = negctrl_data_preproc,
                                                       trt_var = potential_unmeasured_confounder_varname,
                                                       out_var = get_outcome_colname(),
                                                       expected_trt_eff = 0.0)

            if (FALSE == negctrl_outcome_u_verified) {
                warning(glue("!! ERROR: Sim #{sim_number}: Negative control to U1 association is not zero."))
            }
            print(glue("> DEBUG: Sim #{sim_number}: Negative control to U1 association ",
                       "in first negative outcome {negctrl_outcome_u_verified}"))
        } # End checking for zero association between negative control outcomes to U

    } # End verifying first negative control

    ##
    # Negative controls: Creating negative control data -----------------------
    ##

    print("> DEBUG: Preprocessing data to generate OBSERVED negative control")

    trt_data_for_negctrl <- trt_data
    if (TRUE == ignore_extra_confounder_in_negctrl) {
        print(glue("> DEBUG: Sim #{sim_number}: Generating negative control data: removing `{potential_unmeasured_confounder_varname}`"))
        trt_data_for_negctrl <- trt_data_for_negctrl %>%
            select(-.data[[potential_unmeasured_confounder_varname]])
    }

    # Observed data don't have quad term
    trt_data_for_negctrl <- trt_data_for_negctrl %>%
        if_cond_do(condition = (quad_term_sim_scenario & ignore_x1sqr_negctrl_model),
            fun = { function(data) {
                print(glue("> DEBUG: Sim #{sim_number}: ",
                    "Generating negative control data: ",
                    "removing `{get_quadterm_colname()}`"))

                data %>% select(-.data[[get_quadterm_colname()]])
            }
        }) %>% # Observed data don't have X1*X2 term
        if_cond_do(
            condition = (intr_twoconfounders_sim_scenario & ignore_x1x2_in_obsest), fun = { function(data) {
                print(glue("> DEBUG: Sim #{sim_number}: Generating negative control data: ",
                           "removing `{get_x1x2intr_colname()}`"))
                data %>% select(-.data[[get_x1x2intr_colname()]])
        }})

    # Observed data don't have X1*Z term
    trt_data_for_negctrl <- trt_data_for_negctrl %>%
        if_cond_do(
            condition = (intr_conftrt_sim_scenario & ignore_intr_conftrt_in_obsest),
            fun = { function(data) {
                print(glue("> DEBUG: Sim #{sim_number}: ",
                           "Generating negative control data: ",
                           "removing `{get_x1_zintr_colname()}`"))

                       data %>% select(-.data[[get_x1_zintr_colname()]])
                }
            })

    negctrl_outcome_nums <- seq.int(from = 2L, to = num_negctrls + 1L)
    negctrl_data <- negctrl_outcome_nums %>% purrr::map_dfr(.f = create_negctrl_data,
                                                            negctrl_outcomes_mat = outcomes_mat,
                                                            covars_data = trt_data_for_negctrl,
                                                            outcome_var_name = get_outcome_colname())


    if (debug_intrpt) {
        print("> DEBUG: Sim #{sim_number}: -Negative- control data to model `negctrl_data`")
        print(negctrl_data) ;

        print(glue("> DEBUG: Sim #{sim_number}: ",
            "`ignore_extra_confounder_in_negctrl` : {ignore_extra_confounder_in_negctrl}"))
        if (quad_term_sim_scenario) {
            print(glue("> DEBUG: Sim #{sim_number}: ",
                "`ignore_x1sqr_negctrl_model` : {ignore_x1sqr_negctrl_model}"))
        }

        print("> DEBUG: Sim #{sim_number}: First -negative- control data to model")
        print(head(negctrl_data$negctrl_data[[1]], 5))
        browser()
    }
    if (persist_sim_out_data) {
        sim_output[['negctrl_data']] <- negctrl_data
    }

    ##
    # Negative controls: Modelling --------------------------------------------
    ##

    if (num_negctrls < 1) { stop("!! ERROR: Needs to have more than one negative control.") }

    negctrl_outcome_nums <- seq.int(from = 2L, to = num_negctrls + 1L)

    print(glue("> DEBUG: Sim #{sim_number}: Fitting regression to -negative- controls"))

    negctrl_coef_est <- negctrl_data %>% apply(MARGIN = 1, FUN = function(negctrl_data_tbl_row) {

        outcome_num <- negctrl_data_tbl_row[['outcome_id']]
        negctrl_num <- negctrl_data_tbl_row[['negctrl_id']]
        negctrl_data <- negctrl_data_tbl_row[['negctrl_data']]

        print(glue("> DEBUG: Sim #{sim_number}: Fitting IPW GLM to negative control ",
                   "{negctrl_num} (outcome num {outcome_num})"))

        vars_to_exclude_trt_model <- get_trt_colname() %>%
            if_cond_do(condition = (quad_term_sim_scenario & (FALSE == x1sqr_in_trt_model)),
                       fun = { function(vars_to_remove_vec) {
                    vctrs::vec_c(vars_to_remove_vec, get_quadterm_colname())
                }
            }) %>%
            if_cond_do(condition = intr_conftrt_sim_scenario,
                       fun = { function(vars_to_remove_vec) {
                           vctrs::vec_c(vars_to_remove_vec, get_x1_zintr_colname())
                           }
                       })

        negctrl_coef_est <- est_trtcoef_w_sweights(negctrl_data,
                                                   trt_var = get_trt_colname(),
                                                   outcome_var = get_outcome_colname(),
                                                   vars_to_remove_from_trt_model = vars_to_exclude_trt_model,
                                                   vars_to_remove_from_outcome_model = c())

        return(negctrl_coef_est)
    })

    negctrl_coef_est_tbl <- tibble::tibble(
        outcome_id = negctrl_data %>% pull(outcome_id),
        negctrl_id = negctrl_data %>% pull(negctrl_id),
        true_trt_eff = rep(0.0, times = nrow(negctrl_data)),
        est = negctrl_coef_est)

    if (debug_intrpt) {
        print("Estimates using negctrls `negctrl_coef_est_tbl`")
        print(negctrl_coef_est_tbl)
        browser()
    }
    if (TRUE == calibrate_negctrls) {
        sim_output[['negctrl_est']] <- negctrl_coef_est_tbl
    }

    ##
    # Positive controls: Generation -------------------------------------------
    ##

    # Create template of coef table to generate positive controls
    posctrl_coef_prefix = "c"
    posctrl_coef_template <- convert_negctrl_est_to_coef_tbl(
        negctrl_est_tbl = negctrl_coef_est_tbl,
        num_negctrl = num_negctrl,
        trt_covarname = get_trt_colname(),
        confounder_est_prefix = get_confounder_colname_prefix(),
        new_coef_prefix = posctrl_coef_prefix)

    print(glue("> DEBUG: Sim #{sim_number}: Coefficients for outcome model of positive control"))
    if (debug_intrpt) { print(posctrl_coef_template) ; browser() }

    postctrl_data <- negctrl_data %>% purrr::pmap_dfr(.f = ~ {
        outcome_num <- ..1
        negctrl_num <- ..2
        current_negctrl_data <- ..3

        print(glue("> DEBUG: Sim #{sim_number}: Generating positive control using neg ctrl number {negctrl_num} ",
                   "(outcome id {outcome_num})"))

        # Select appropriate row of coefficient template
        posctrl_coef_template_tbl_row <- ..4 %>%
            filter(.data[['outcome_id']] == outcome_num) %>%
            select(-.data[['outcome_id']])
        coef_vec_names <- names(posctrl_coef_template_tbl_row)

        posctrl_coef_template_vec <- posctrl_coef_template_tbl_row %>% slice(1) %>% as.numeric()
        names(posctrl_coef_template_vec) <- coef_vec_names

        coef_prefix <- ..5
        outcome_varname <- ..6
        trt_varname <- ..7

        print(glue("> DEBUG: Sim #{sim_number} NegCtrl #{negctrl_num}: Adjusting positive control ",
                   "true effect size by subtracting treatment estimate of negative control"))
        negctrl_trt_eff_est <- negctrl_coef_est_tbl$est[[negctrl_num]] %>%
            filter(term == get_trt_colname()) %>% pull(estimate)
        adj_posctrl_log_or_targets <- posctrl_log_or_targets - negctrl_trt_eff_est

        current_posctrl_data <- adj_posctrl_log_or_targets %>% purrr::map_dfr(
            .f = gen_posctrl_data_v2,
            data = current_negctrl_data,
            coef_template_vec = posctrl_coef_template_vec,
            coef_prefix = coef_prefix,
            trt_varname = trt_varname, outcome_varname = outcome_varname)

        # Adds additional information
        current_posctrl_data <- current_posctrl_data %>%
            mutate(outcome_id = rep(outcome_num, times = nrow(current_posctrl_data))) %>%
            mutate(negctrl_id = rep(negctrl_num, times = nrow(current_posctrl_data))) %>%
            relocate(negctrl_id) %>% relocate(outcome_id)

    }, posctrl_coef_template_tbl = posctrl_coef_template,
       posctrl_coef_prefix = posctrl_coef_prefix,
       outcome_varname = get_outcome_colname(), trt_varname = get_trt_colname() )

    if (debug_intrpt) { browser() }
    postctrl_data <- postctrl_data %>%
        mutate(true_trt_eff = rep(posctrl_log_or_targets, times = num_negctrls))

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: Positive control data `postctrl_data`"))
        print(postctrl_data)
        browser()
    }

    ##
    # Positive controls: Modelling --------------------------------------------
    ##

    print(glue("> DEBUG: Sim #{sim_number}: Modelling positive controls"))
    vars_to_exclude_trt_model <- get_trt_colname() %>% if_cond_do(
        condition = (quad_term_sim_scenario & (FALSE == x1sqr_in_trt_model)), fun = {
            function(vars_to_remove_vec) {
                    vctrs::vec_c(vars_to_remove_vec, get_quadterm_colname())
            }
        })
    vars_to_exclude_trt_model <- c(vars_to_exclude_trt_model, "true_trt_eff")
    vars_to_exclude_outcome_model <- c("true_trt_eff")

    if (debug_intrpt) { browser() }

    postctrl_est_tbl <- postctrl_data %>% purrr::pmap_dfr(.f = ~{
        outcome_id <- ..1
        negctrl_id = ..2
        posctrl_data <- ..4
        est_coef_posctrl(outcome_id, negctrl_id, posctrl_data,
                         outcome_varname = get_outcome_colname(),
                         trt_varname = get_trt_colname(),
                         vars_to_remove_from_trt_model = vars_to_exclude_trt_model,
                         vars_to_remove_from_outcome_model = vars_to_exclude_outcome_model,
                         use_ipw = TRUE)
    })
    postctrl_est_tbl <- postctrl_est_tbl %>%
        mutate(true_trt_eff = postctrl_data %>% pull(true_trt_eff)) %>%
        relocate(true_trt_eff, .before = est)
    postctrl_est_tbl <- postctrl_est_tbl %>%
        mutate(adj_true_trt_eff = postctrl_data %>% pull(trt_eff))

    if (debug_intrpt) { print(postctrl_est_tbl) ; browser() }

    ##
    # Empirical calibration: Collect data to build the systematic error model -----------------
    ##

    # Collate estimates from negative and positive controls
    #
    # `posctrl_has_adj_true_eff` flag adjusts for
    # Systematic error model with ADJUSTED POSITIVE controls

    print(glue("> DEBUG: Sim #{sim_number}: Collating estimates from ALL controls"))
    all_ctrl_est_tbl <- collate_est_from_ctrls(negctrl_est_tbl = negctrl_coef_est_tbl,
                                               posctrl_est_tbl = postctrl_est_tbl,
                                               trt_varname = get_trt_colname(),
                                               posctrl_has_adj_true_eff = TRUE,
                                               ctrl_tbl_colnames = list(
                                                   outcome_id_colname = 'outcome_id',
                                                   negctrl_id_colname = 'negctrl_id',
                                                   est_colname = 'est',
                                                   true_trt_eff_colname = 'true_trt_eff',
                                                   adj_true_trt_eff_colname = 'adj_true_trt_eff'))

    all_ctrl_est_tbl <- all_ctrl_est_tbl %>%
        mutate(sim_id = rep(sim_number, times = nrow(all_ctrl_est_tbl))) %>%
        relocate(sim_id)

    all_ctrl_est_tbl <- all_ctrl_est_tbl %>% arrange(sim_id, outcome_id, true_eff)
    if (debug_intrpt) { print(all_ctrl_est_tbl) ; browser() }
    sim_output[['all_ctrl_est_tbl']] <- all_ctrl_est_tbl


    ## NEGATIVE controls only data ------
    negctrl_calib_est <- all_ctrl_est_tbl %>% filter(true_eff == 0)

    ##
    # Empirical calibration: Build systematic error model ----------------------
    ##

    if (debug_intrpt & (sim_num_to_plot == sim_number)) {
        sys_err_model_plot <- EmpiricalCalibration::plotErrorModel(
            logRr = all_ctrl_est_tbl[['actual_eff_est']],
            seLogRr = all_ctrl_est_tbl[['actual_eff_est_stderr']],
            trueLogRr = all_ctrl_est_tbl[['true_eff']])
        if (debug_intrpt) { print(sys_err_model_plot) ; browser() }
    }

    if (debug_intrpt & (sim_num_to_plot == sim_number)) {
        ctrl_est_forest_plot <- plotTrueAndObserved(all_ctrl_est_tbl[['actual_eff_est']],
                                                    all_ctrl_est_tbl[['actual_eff_est_stderr']],
                                                    all_ctrl_est_tbl[['true_eff']])
        ctrl_est_forest_plot <- ctrl_est_forest_plot + xlab("True effect size") +
            ggtitle("Forest plot of effect sizes from controls",
                    subtitle = glue("Sim #{sim_number}. Orange: Statistically significant estimate"))

        if (debug_intrpt) { print(ctrl_est_forest_plot) ; browser() }
    }

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: Sys Err Model Legacy flag : {sys_err_model_legacy_flag}"))
        browser()
    }

    print(glue("> DEBUG: Sim #{sim_number}: Creating the systematic error model"))
    sys_err_model <- EmpiricalCalibration::fitSystematicErrorModel(
        logRr = all_ctrl_est_tbl[['actual_eff_est']],
        seLogRr = all_ctrl_est_tbl[['actual_eff_est_stderr']],
        trueLogRr = all_ctrl_est_tbl[['true_eff']],
        legacy = sys_err_model_legacy_flag)

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: Sys Err Model"))
        print(sys_err_model) ; browser()
    }
    sim_output[['sys_err_allctrls']] <- sys_err_model

    ## Systematic error model by converting the NULL model. -------------------
    null_err_model <- EmpiricalCalibration::fitNull(
        logRr = negctrl_calib_est[['actual_eff_est']],
        seLogRr = negctrl_calib_est[['actual_eff_est_stderr']]
    )
    sys_err_model_from_null <- EmpiricalCalibration::convertNullToErrorModel(
        null_err_model, meanSlope = 1, sdSlope = 0) # Defaults

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: Sys Err Model using converted NULL model"))
        print(null_err_model) ; browser()
    }
    sim_output[['sys_err_from_null']] <- sys_err_model_from_null

    ## Systematic error model using NEGATIVE CONTROLS -------------------------

    #sys_err_model_negctrl_only <- EmpiricalCalibration::fitSystematicErrorModel(
    #    logRr     = negctrl_calib_est[['actual_eff_est']],
    #    seLogRr   = negctrl_calib_est[['actual_eff_est_stderr']],
    #    trueLogRr = negctrl_calib_est[['true_eff']],
    #    legacy = sys_err_model_legacy_flag)

    # Systematic error model negative control is the same as NULL model.
    sys_err_model_negctrl_only <- sys_err_model_from_null

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: Sys Err Model using NEGATIVE CONTROLS only"))
        print(sys_err_model_negctrl_only) ; browser()
    }
    sim_output[['sys_err_model_negctrl_only']] <- sys_err_model_negctrl_only

    ## Systematic error model using ADJUSTED POSTITIVE CONTROLS ---------------
    sys_err_model_adj_true_eff <- EmpiricalCalibration::fitSystematicErrorModel(
        logRr = all_ctrl_est_tbl[['actual_eff_est']],
        seLogRr = all_ctrl_est_tbl[['actual_eff_est_stderr']],
        trueLogRr = all_ctrl_est_tbl[['adj_true_eff']],
        legacy = sys_err_model_legacy_flag)

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: Sys Err Model using adjusted TRUE effects of POSITIVE CONTROLS"))
        print(sys_err_model_adj_true_eff) ; browser()
    }
    sim_output[['sys_err_model_adj_true_eff']] <- sys_err_model_adj_true_eff

    ##
    # Outcome of interest: Modelling ------------------------------------------
    ##

    print(glue("> DEBUG: Sim #{sim_number}: Extracting data for *outcome of interest*"))
    tmp_data <-
        create_negctrl_data(outcome_id = 1L,
                            negctrl_outcomes_mat = outcomes_mat,
                            covars_data = trt_data,
                            outcome_var_name = get_outcome_colname())
    outcome_of_interest_data <- tmp_data$negctrl_data[[1]]

    sim_output[['primary_outcome_data']] <- outcome_of_interest_data

    if (debug_intrpt) {
        print("> DEBUG: Sim #{sim_number}: Data for primary outcome `outcome_of_interest_data`")
        print(head(outcome_of_interest_data, 5)) ; browser()
    }

    print(glue("> DEBUG: Sim #{sim_number}: Estimating treatment effect for *outcome of interest*"))

    vars_to_remove_from_outcome_model <- c(get_outcome_colname()) %>%
        if_cond_do(condition = ignore_extra_confounder_in_study,
            fun = { function(vars_to_rm_vec) {
                print(glue("> DEBUG: Sim #{sim_number}: Modelling ",
                    "outcome of interest: ignoring U"))
                vctrs::vec_c(vars_to_rm_vec, potential_unmeasured_confounder_varname)
            }
        }) %>%
        if_cond_do(condition = (quad_term_sim_scenario & ignore_x1sqr_study_model),
            fun = { function(vars_to_rm_vec) {
                print(glue("> DEBUG: Sim #{sim_number}: Modelling ",
                    "outcome of interest: removing X1Sqr"))
                vctrs::vec_c(vars_to_rm_vec, get_quadterm_colname())
            }
        })

    vars_to_remove_from_trt_model <- c(get_trt_colname()) %>%
            if_cond_do(condition = ignore_extra_confounder_in_study,
            fun = { function(vars_to_rm_vec) {
                print(glue("> DEBUG: Sim #{sim_number}: Modelling ",
                    "TREATMENT in primary outcome: removing U"))
                vctrs::vec_c(vars_to_rm_vec, potential_unmeasured_confounder_varname)
            }
        }) %>%
        if_cond_do(condition = (quad_term_sim_scenario & (FALSE == x1sqr_in_trt_model)),
            fun = { function(vars_to_rm_vec) {
                print(glue("> DEBUG: Sim #{sim_number}: Modelling ",
                    "TREATMENT in primary outcome: removing X1Sqr"))
                vctrs::vec_c(vars_to_rm_vec, get_quadterm_colname())
            }
        })

    outcome_of_intest_est <- est_trtcoef_w_sweights(outcome_of_interest_data,
                           trt_var = get_trt_colname(), outcome_var = get_outcome_colname(),
                           vars_to_remove_from_trt_model = vars_to_remove_from_trt_model,
                           vars_to_remove_from_outcome_model = vars_to_remove_from_outcome_model)

    if (debug_intrpt) {
        print(outcome_of_intest_est)
        print(glue(">> DEBUG: Sim #{sim_number}: True treatment effect for primary outcome ",
                   "{signif(trt_eff, 3)}"))
        browser()
    }

    ##
    # Empirical calibration: Calibrate outcome of interest ----------------------------------------
    ##

    print(glue("> DEBUG: Sim #{sim_number}: ",
               "Calibrating CI of outcome of interest `calibrated_est`"))
    calibrated_est <- EmpiricalCalibration::calibrateConfidenceInterval(
        outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(estimate),
        outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(std.error),
        sys_err_model)

    if (debug_intrpt) { print(calibrated_est) ; browser() }

    ## Calibrate with NEGATIVE CONTROLS only -----
    print(glue("> DEBUG: Sim #{sim_number}: Calibrating using NEGATIVE CONTROLS"))
    calibrated_est_negctrl_only <- EmpiricalCalibration::calibrateConfidenceInterval(
        outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(estimate),
        outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(std.error),
        sys_err_model_negctrl_only)

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: ",
                   "Calibrated CI of outcome of interest using neg controls only `calibrated_est_negctrl_only`"))
        print(calibrated_est_negctrl_only)
        browser()
    }

    ## Calibrate with converted NULL model ----
    print(glue("> DEBUG: Sim #{sim_number}: Calibrating using converted NULL model"))
    calibrated_est_null_model <- EmpiricalCalibration::calibrateConfidenceInterval(
        outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(estimate),
        outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(std.error),
        sys_err_model_from_null)

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: ",
                   "Calibrated CI of outcome of interest using converted null model, `calibrated_est_null_model`"))
        print(calibrated_est_null_model)
        browser()
    }

    ## Calibrate with ADJUSTED TRUE EFFECT OF POSITIVE CONTROLS -----
    print(glue("> DEBUG: Sim #{sim_number}: Calibrating using ",
        "adjusted true effect of POSTIVIE CONTROLS"))
    calibrated_est_adjtrueff <- EmpiricalCalibration::calibrateConfidenceInterval(
        outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(estimate),
        outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(std.error),
        sys_err_model_adj_true_eff)

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: ",
                   "Calibrated CI of outcome of interest using ",
                   "adjusted truth of neg controls only `calibrated_est_adjtrueff`"))
        print(calibrated_est_adjtrueff)
        browser()
    }

    ##
    # Calculates biases -------------------------------------------------------
    ##

    print(glue("> DEBUG: Sim #{sim_number}: True treatment effect {signif(true_trt_eff_persim, 3)}"))

    # Obtains the names of systematic error models
    sys_err_types <- get_syserr_model_calibration_type()

    bias_calc_tbl <- tibble(
        calitype = rep(sys_err_types[['nomod']], times = 2),
        restype = c("uncalibrated", "calibrated"),
        est = c(outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(estimate),
                calibrated_est[['logRr']]),
        est_stderr = c(outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(std.error),
                       calibrated_est[['seLogRr']]),
        true_eff = rep(true_trt_eff_persim, times = 2))
    bias_calc_tbl <- bias_calc_tbl %>% mutate(est_bias = true_trt_eff_persim - est)

    sim_output[['bias_calc']] <- bias_calc_tbl

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: Bias calculation `bias_calc_tbl`"))
        print(bias_calc_tbl)
        browser()
    }

    if (sum(is.na(bias_calc_tbl %>% pull(est_stderr))) > 0) {
        warning(glue("!! WARNING: Sim #{sim_number}: Standard Error of calibrated result is NA"))

        print(all_ctrl_est_tbl)
        sys_err_model_plot <- EmpiricalCalibration::plotErrorModel(
            all_ctrl_est_tbl[['actual_eff_est']],
            all_ctrl_est_tbl[['actual_eff_est_stderr']],
            all_ctrl_est_tbl[['true_eff']])
        print(sys_err_model_plot)
        browser()
    }

    ## Bias calculation: error model constructed with neg ctrl only --------
    bias_negctrl_only_sys_err_tbl <- tibble(
        calitype = rep(sys_err_types[['negctrlonlysyserr']], times = 2),
        restype = c("uncalibrated", "calibrated"),
        est = c(outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(estimate),
                calibrated_est_negctrl_only[['logRr']]),
        est_stderr = c(outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(std.error),
                       calibrated_est_negctrl_only[['seLogRr']]),
        true_eff = rep(true_trt_eff_persim, times = 2))

    bias_negctrl_only_sys_err_tbl <- bias_negctrl_only_sys_err_tbl %>% mutate(est_bias = true_trt_eff_persim - est)
    sim_output[['bias_calc_negctrlonlysyserr']] <- bias_negctrl_only_sys_err_tbl

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: Bias calculation using ",
                   " negctrl only `bias_negctrl_only_sys_err_tbl`"))
        print(bias_negctrl_only_sys_err_tbl)
        browser()
    }

    ## Bias calculation: error model constructed with converted NULL only --------
    bias_null_sys_err_tbl <- tibble(
        calitype = rep(sys_err_types[['nullsyserr']], times = 2),
        restype = c("uncalibrated", "calibrated"),
        est = c(outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(estimate),
                calibrated_est_null_model[['logRr']]),
        est_stderr = c(outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(std.error),
                       calibrated_est_null_model[['seLogRr']]),
        true_eff = rep(true_trt_eff_persim, times = 2))

    bias_null_sys_err_tbl <- bias_null_sys_err_tbl %>% mutate(est_bias = true_trt_eff_persim - est)
    sim_output[['bias_calc_null_syserr']] <- bias_null_sys_err_tbl

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: Bias calculation using ",
                   " null model. `bias_null_sys_err_tbl` :"))
        print(bias_null_sys_err_tbl)
        browser()
    }

    ## Bias calculation: error model constructed with adjusted truth of positive controls -----
    bias_adjtrueff_sys_err_tbl <- tibble(
        calitype = rep(sys_err_types[['adjtruthposctrlsyserr']], times = 2),
        restype = c("uncalibrated", "calibrated"),
        est = c(outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(estimate),
                calibrated_est_adjtrueff[['logRr']]),
        est_stderr = c(outcome_of_intest_est %>% filter(term == get_trt_colname()) %>% pull(std.error),
                       calibrated_est_adjtrueff[['seLogRr']]),
        true_eff = rep(true_trt_eff_persim, times = 2))
    bias_adjtrueff_sys_err_tbl <- bias_adjtrueff_sys_err_tbl %>% mutate(est_bias = true_trt_eff_persim - est)
    sim_output[['bias_calc_adjtruthposctrlsyserr']] <- bias_adjtrueff_sys_err_tbl

    if (debug_intrpt) {
        print(glue("> DEBUG: Sim #{sim_number}: Bias calculation using ",
                   " adjusted truth in positive controls `bias_adjtrueff_sys_err_tbl`"))
        print(bias_adjtrueff_sys_err_tbl)
        browser()
    }

    ##
    # Calibrates the CONTROLS Suchard, et al. Lancet --------------------------
    ##

    ## Calibrating the negative controls ---------------------------------------
    if (TRUE == calibrate_negctrls) {

        print("> DEBUG: Calibrating negative controls")
        if (debug_intrpt) { browser() }

        ### Calibrate negative controls using systematic error that contains all controls ----

        ctrl_data_to_calibrate <- all_ctrl_est_tbl %>% filter(true_eff == 0.0)

        calibrated_negctrl_data <- calibrate_control(est_from_all_ctrls_tbl = ctrl_data_to_calibrate,
                                                     sys_err_model = sys_err_model,
                                                     all_ctrls_tbl_coldef = list(true_eff_colname = "true_eff",
                                                                                       est_colname = "actual_eff_est",
                                                                                       se_est_colname = "actual_eff_est_stderr"))

        calibrated_negctrl_est_tbl_v2 <- merge_calib_uncalib_ctrl_data(
            ctrl_data_to_calibrate, calibrated_negctrl_data,
            cali_sys_err_model_type = get_syserr_model_calibration_type()[['nomod']])

        if (debug_intrpt) { browser() }

        #sim_output[['calibrated_negctrl_est_old']] <- calibrated_negctrl_est_tbl
        sim_output[['calibrated_negctrl_est']] <- calibrated_negctrl_est_tbl_v2

        ### Calibrate negative controls using systematic error that contains only negative controls ----
        print("> DEBUG: Calibrating negative controls using sys err of neg controls only")
        calibrated_negctrlonly_negctrl_data <- calibrate_control(
                                                     est_from_all_ctrls_tbl = ctrl_data_to_calibrate,
                                                     sys_err_model = sys_err_model_negctrl_only,
                                                     all_ctrls_tbl_coldef = list(true_eff_colname = "true_eff",
                                                                                 est_colname = "actual_eff_est",
                                                                                 se_est_colname = "actual_eff_est_stderr"))
        calibrated_negctrl_est_tbl_v3 <- merge_calib_uncalib_ctrl_data(
            ctrl_data_to_calibrate, calibrated_negctrlonly_negctrl_data,
            cali_sys_err_model_type = get_syserr_model_calibration_type()[['negctrlonlysyserr']])

        if (debug_intrpt) { browser() }
        sim_output[['calibrated_negctrl_est_negctrlonlysyserr']] <- calibrated_negctrl_est_tbl_v3

        ### Calibrate negative controls using converted NULL systematic error ----
        print("> DEBUG: Calibrating negative controls using convereted null sys err model")
        calibrated_nullmodel_negctrl_data <- calibrate_control(est_from_all_ctrls_tbl = ctrl_data_to_calibrate,
                                                               sys_err_model = sys_err_model_from_null,
                                                               all_ctrls_tbl_coldef = list(true_eff_colname = "true_eff",
                                                                                             est_colname = "actual_eff_est",
                                                                                             se_est_colname = "actual_eff_est_stderr"))
        calibrated_negctrl_est_tbl_v5 <- merge_calib_uncalib_ctrl_data(
            ctrl_data_to_calibrate, calibrated_nullmodel_negctrl_data,
            cali_sys_err_model_type = get_syserr_model_calibration_type()[['nullsyserr']])
        if (debug_intrpt) { browser() }
        sim_output[['calibrated_negctrl_est_nullsyserr']] <- calibrated_negctrl_est_tbl_v5


        ### Calibrate negative controls using systematic error that contains only ADJUSTED positive controls ----
        print("> DEBUG: Calibrating negative controls using sys err of ADJUSTED truth positive controls")
        calibrated_w_adjposctrl_negctrl_data <- calibrate_control(est_from_all_ctrls_tbl = ctrl_data_to_calibrate,
                                                                 sys_err_model = sys_err_model_adj_true_eff,
                                                                 all_ctrls_tbl_coldef = list(true_eff_colname = "adj_true_eff",
                                                                                             est_colname = "actual_eff_est",
                                                                                             se_est_colname = "actual_eff_est_stderr"))
        calibrated_negctrl_est_tbl_v4 <- merge_calib_uncalib_ctrl_data(
            ctrl_data_to_calibrate, calibrated_w_adjposctrl_negctrl_data,
            cali_sys_err_model_type = get_syserr_model_calibration_type()[['adjtruthposctrlsyserr']])

        if (debug_intrpt) { browser() }
        sim_output[['calibrated_negctrl_est_adjtruesyserr']] <- calibrated_negctrl_est_tbl_v4


    } # End calibrating negative controls

    ## Calibrating the positive controls ---------------------------------------
    if (TRUE == calibrate_posctrls) {
        print("> DEBUG: Calibrating positive controls")
        if (debug_intrpt) { browser() }

        posctrl_to_calibrate <- all_ctrl_est_tbl %>% filter(true_eff != 0.0)

        print("> DEBUG: Calibrating positive controls using all controls")
        calibrated_posctrl_data_df <- posctrl_to_calibrate %>%
            calib_posctrl(err_model = sys_err_model,
                          use_adj_true_colname = FALSE)
        sim_output[['calibrated_posctrl_est']] <- calibrated_posctrl_data_df %>%
            mutate(calitype = rep(get_syserr_model_calibration_type()[['nomod']],
                                  times = nrow(calibrated_posctrl_data_df)))

        if (debug_intrpt) { browser() }

        print("> DEBUG: Calibrating positive controls using negative controls only")
        calibrated_posctrl_w_negctrl_only_data_df <- posctrl_to_calibrate %>%
            calib_posctrl(err_model = sys_err_model_negctrl_only,
                          use_adj_true_colname = FALSE)
        sim_output[['calibrated_posctrl_est_negctrlonlysyserr']] <- calibrated_posctrl_w_negctrl_only_data_df %>%
            mutate(calitype =  rep(get_syserr_model_calibration_type()[['negctrlonlysyserr']],
                                  times = nrow(calibrated_posctrl_w_negctrl_only_data_df)))

        if (debug_intrpt) { browser() }

        print("> DEBUG: Calibrating positive controls using NULL model")
        calibrated_posctrl_w_nullmodel_data_df <- posctrl_to_calibrate %>%
            calib_posctrl(err_model = sys_err_model_from_null,
                          use_adj_true_colname = FALSE)
        sim_output[['calibrated_posctrl_est_nullsyserr']] <- calibrated_posctrl_w_nullmodel_data_df %>%
            mutate(calitype = rep(get_syserr_model_calibration_type()[['nullsyserr']],
                                  times = nrow(calibrated_posctrl_w_negctrl_only_data_df)))

        if (debug_intrpt) { browser() }

        print("> DEBUG: Calibrating positive controls using sys err of ADJUSTED truth positive controls")
        calibrated_posctrl_w_adjtrue_df <- posctrl_to_calibrate %>%
            calib_posctrl(err_model = sys_err_model_adj_true_eff,
                          use_adj_true_colname = TRUE)
        sim_output[['calibrated_posctrl_est_adjtruesyserr']] <- calibrated_posctrl_w_adjtrue_df %>%
            mutate(calitype = rep(get_syserr_model_calibration_type()[['adjtruthposctrlsyserr']],
                                  times = nrow(calibrated_posctrl_w_negctrl_only_data_df)))

        if (debug_intrpt) { browser() }

    } # End (TRUE == calibrate_posctrls)

    return(sim_output)

} # End of simulation

if (n_sim > parallelise_sim_threshold) {
    parallel::stopCluster(cluster)
    registerDoSEQ()
}

print("> Finished simulation")  # Simulation completes -------------------------
sim_complete_time <- lubridate::now()

sim_interval <- lubridate::interval(sim_start_time, end = sim_complete_time)
sim_duration_min <- lubridate::int_length(sim_interval) / 60

print(sim_complete_time)

##
# Parsing results from simulations --------------------------------------------
##

sim_complete_time_formatted <- sim_complete_time %>% format("%Y%m%d_%H%M")

print("> DEBUG: Creating experiment settings")

same_all_measured_eff_to_all_outcomes_all_sim_val <- NA
same_measured_eff_to_all_outcomes_per_sim_iter_val <- NA
common_coef_all_confounders_all_outcomes_val <- NA


if (FALSE == use_outcome_model_coef_template_all_sims) {
    same_all_measured_eff_to_all_outcomes_all_sim_vals <- extract_tbl_from_sim_result(sim_res = sim_out, name = "same_all_measured_eff_to_all_outcomes_all_sim")
    same_all_measured_eff_to_all_outcomes_all_sim_val <- same_all_measured_eff_to_all_outcomes_all_sim_vals %>% pull(res) %>% unlist() %>% unique()

    same_measured_eff_to_all_outcomes_per_sim_iter_vals <- extract_tbl_from_sim_result(sim_res = sim_out, name = "same_measured_eff_to_all_outcomes_per_sim_iter")
    same_measured_eff_to_all_outcomes_per_sim_iter_val <- same_measured_eff_to_all_outcomes_per_sim_iter_vals %>% pull(res) %>% unlist() %>% unique()

    common_coef_all_confounders_all_outcomes_vals <- extract_tbl_from_sim_result(sim_res = sim_out, name = "common_coef_all_confounders_all_outcomes")
    common_coef_all_confounders_all_outcomes_val <- common_coef_all_confounders_all_outcomes_vals %>% pull(res) %>% unlist() %>% unique()
}

num_outcomes <- num_negctrls + 1L # Plus one for outcome of interest

sim_settings <- tribble(
    ~flags                              , ~persim_allsim_flags,           ~value,
    'all_sim_use_trt_model_template'    ,            'AllSims', all_sim_use_trt_model_template,
    'gen_extra_confounder_in_study'   ,            'AllSims', gen_extra_confounder_in_study,
    'gen_extra_confounder_in_negctrl' ,            'AllSims', gen_extra_confounder_in_negctrl,
    'ignore_extra_confounder_in_study',             'PerSim', ignore_extra_confounder_in_study,
    'ignore_extra_confounder_in_negctrl',           'PerSim', ignore_extra_confounder_in_negctrl,
    'use_outcome_model_coef_template_all_sims',      'AllSims', use_outcome_model_coef_template_all_sims,
    'use_same_extra_term_coef_vec_all_sim',          'AllSims', use_same_extra_term_coef_vec_all_sim,
    'outcome_model_coef_confounder_same_eff_all_outcomes_all_sims', 'AllSims', outcome_model_coef_confounder_same_eff_all_outcomes_all_sims,
    'confounder_to_all_outcomes_all_sims_eff',       'AllSims', confounder_to_all_outcomes_all_sims_eff,
    'equi_confounding',                              'AllSims', equi_confounding,
    'same_all_measured_eff_to_all_outcomes_all_sim', 'AllSims', same_all_measured_eff_to_all_outcomes_all_sim_val,
    'same_measured_eff_to_all_outcomes_per_sim_iter','AllSims', same_measured_eff_to_all_outcomes_per_sim_iter_val,
    # 'common_coef_all_confounders_all_outcomes_val',  'AllSims', common_coef_all_confounders_all_outcomes_val, # List?
    'n_sims',                                         'PerSim', n_sim,
    'n_obs_persim',                                   'PerSim', n_obs,
    'n_negctrls',                                     'PerSim', num_negctrls,
    'num_not_trt_pri_outcome_mconfounder',           'AllSims', num_notpriout_mconf,
    'sim_duration_mins',                             'AllSims', sim_duration_min,
    'calibrate_negctrls',                            'AllSims', calibrate_negctrls,
    'calibrate_posctrls',                            'AllSims', calibrate_posctrls,
    'sys_err_model_legacy' ,                         'AllSims', sys_err_model_legacy_flag
)
print(sim_settings)


# Persisting data ----

# Check and create required directories
if (persist_sim_out_data) {
  if (FALSE == dir.exists(here("data", "sim-binaryoutcome"))) {
      dir.create(here("data", "sim-binaryoutcome"))
  }  
}
if (persist_plots) {
  if (FALSE == dir.exists(here("graphs"))) {
      dir.create(here("graphs"))
  }
}

if (persist_sim_out_data) {
    sim_settings_file <- here::here("data", "sim-binaryoutcome",
                                      glue("simres-{sim_complete_time_formatted}-settings-{sim_name}.csv"))
    print( glue("> DEBUG: Persisting sim settings to `{sim_settings_file}`"))
    readr::write_csv(sim_settings,sim_settings_file)
}

## Persisting randomised measurement error parameters ----
if (measurement_err_confounder_scenario & rand_merror_gaussian_params & persist_sim_out_data) {
    print("> DEBUG: Extracting randomised measurement error parameters")

    merror_gaussian_params_raw <- extract_tbl_from_sim_result(sim_res = sim_out, name = "merror_gaussian_params")
    merror_gaussian_params <- merror_gaussian_params_raw %>% tidyr::unnest(res)

    merror_param_file <- here::here("data", "sim-binaryoutcome",
                                      glue("simres_raw-{sim_complete_time_formatted}-merror_gaussparams-{sim_name}.csv"))
    print( glue("> DEBUG: Persisting measurement error parameters to `{merror_param_file}`"))
    readr::write_csv(merror_gaussian_params, merror_param_file)
}


## Persisting bias calculations -----
print("> DEBUG: Extracting bias calculations from simulations")

bias_calc_from_ctrls_raw <- extract_tbl_from_sim_result(sim_res = sim_out, name = "bias_calc")
bias_calc_data_allctrls <- bias_calc_from_ctrls_raw %>% tidyr::unnest(res)

if (persist_sim_out_data) {
    bias_calc_data_file <- here::here("data", "sim-binaryoutcome",
                                      glue("simres_raw-{sim_complete_time_formatted}-bias_raw-{sim_name}.csv"))
    print( glue("> DEBUG: Persisting ALL bias calculation data to `{bias_calc_data_file}`"))
    readr::write_csv(bias_calc_data_allctrls, bias_calc_data_file)
}

### Persisting bias calculation from using negative controls only -----
bias_calc_negctrl_only <- extract_tbl_from_sim_result(sim_res = sim_out,
    name = "bias_calc_negctrlonlysyserr")
bias_calc_data_negctrl_only <- bias_calc_negctrl_only %>% tidyr::unnest(res)

if (persist_sim_out_data) {
    bias_calc_negctrl_only_data_file <- here::here("data", "sim-binaryoutcome",
        glue("simres_raw-{sim_complete_time_formatted}-bias_negctrlonly-{sim_name}.csv"))
    print( glue("> DEBUG: Persisting bias calculation using negative control ",
        "data to `{bias_calc_negctrl_only_data_file}`"))

    readr::write_csv(bias_calc_data_negctrl_only, bias_calc_negctrl_only_data_file)
}


### Combining all the bias estimates -----
calib_result_from_sim <- rbind(bias_calc_data_allctrls,
                            bias_calc_data_negctrl_only)
                            # bias_calc_data_adj_posctrl)

if (persist_sim_out_data) {
    combined_bias_calc_data_file <- here::here("data", "sim-binaryoutcome",
                                      glue("simres-{sim_complete_time_formatted}-bias_combined_all-{sim_name}.csv"))
    print( glue("> DEBUG: Persisting combined bias calculation data to `{combined_bias_calc_data_file}`"))
    readr::write_csv(calib_result_from_sim, combined_bias_calc_data_file)
}

# Removing calibrated estimates with NAs from bias calculation from all ctrls estimate ----
print("> DEBUG: Removing simulations with NA calibrated standard error")

bias_calc_from_ctrls <- bias_calc_from_ctrls_raw %>% tidyr::unnest(c(sim_num, res))

# Removes studies with NA standard errors
no_na_data_obj <- bias_calc_from_ctrls %>% filter_sims_without_na(na_col = "est_stderr", group_var = "sim_num")
num_na <- no_na_data_obj[['na_groups']] %>% length()
# Removes studies with NA estimates errors
no_na_data_obj <- no_na_data_obj %>% .[['data_without_na']] %>% filter_sims_without_na(na_col = "est", group_var = "sim_num")
num_na <- num_na + no_na_data_obj[['na_groups']] %>% length()
# Repack
bias_calc_from_ctrls <- no_na_data_obj[['data_without_na']] %>%
    tidyr::nest(res = c(calitype, restype, est, est_stderr, true_eff, est_bias))
names(bias_calc_from_ctrls) <- c("sim_num", "res")

#browser("Testing what happens when NAs are not removed")
#bias_calc_from_ctrls <- bias_calc_from_ctrls_raw

if (persist_sim_out_data) {
    bias_calc_data_file <- here::here("data", "sim-binaryoutcome",
                                      glue("simres-{sim_complete_time_formatted}-bias-{sim_name}.csv"))
    print( glue("> DEBUG: Persisting bias calculation data to `{bias_calc_data_file}`"))
    readr::write_csv(bias_calc_from_ctrls %>% tidyr::unnest(res),
                     bias_calc_data_file)
}

print("> DEBUG: Extracting reg coef estimations of controls from simulations; `est_from_ctrls`")
est_from_ctrls <- extract_tbl_from_sim_result(sim_res = sim_out, name = "all_ctrl_est_tbl")
print(est_from_ctrls)

if (persist_sim_out_data) {
    est_data_file <- here::here("data", "sim-binaryoutcome",
                                      glue("simres-{sim_complete_time_formatted}-estimates-{sim_name}.csv"))
    print( glue("> DEBUG: Persisting estimates data to `{est_data_file}`"))
    readr::write_csv(est_from_ctrls %>% tidyr::unnest(res),
                     est_data_file)
}

if (persist_sim_out_data) {
    trt_model_coefs <- extract_tbl_from_sim_result(sim_res = sim_out, name = "trt_model_coef")
    trt_model_coefs_file <- here::here("data", "sim-binaryoutcome",
                                       glue("simres-{sim_complete_time_formatted}-trt_model_coefs-{sim_name}.csv"))
    print( glue("> DEBUG: Persiting reg coefs of TREATMENT models to `{trt_model_coefs_file}`"))
    readr::write_csv(trt_model_coefs %>% tidyr::unnest(res), trt_model_coefs_file)


    outcome_model_coefs <- extract_tbl_from_sim_result(sim_res = sim_out, name = "outcome_model_coef")
    outcome_model_coefs_file <- here::here("data", "sim-binaryoutcome",
                                glue("simres-{sim_complete_time_formatted}-outcome_model_coefs-{sim_name}.csv"))
    print( glue("> DEBUG: Persiting reg coefs of OUTCOME models to `{outcome_model_coefs_file}`"))
    readr::write_csv(outcome_model_coefs %>% tidyr::unnest(res), outcome_model_coefs_file)

}

##
# Plotting --------------------------------------------------------------------
##

# Setting colour palette
palette("Okabe-Ito")
# Functional to save plots in A4 paper size.
ggsave_a4 <- partial(ggsave, units = "mm",  width = 210, height = 297)
ggsave_a4r <- partial(ggsave, units = "mm",  width = 297, height = 210)

common_plot_subtitle <- glue(sim_name , "\n", "{num_na} simulations removed due to NA in calibrated SE")

## Box plot of biases ----

bias_calc_boxplot <- plot_bias_calc_boxplot(bias_calc_sims = bias_calc_from_ctrls,
                                 plot_title = "Boxplot of biases",
                                 sim_desc = common_plot_subtitle,
                                 axis_break_gap = 0.025)

bias_calc_boxplot <- zoom_in_bias_calc_boxplot(bias_boxplot = bias_calc_boxplot,
                                               bias_calc_sims = bias_calc_from_ctrls,
                                               quantile_prob_limits = c(0.0, 0.95),
                                               bias_colname = "est_bias",
                                               plot_title = "Boxplot of biases zoomed 0.0-0.95 quantiles",
                                               sim_desc = common_plot_subtitle)
print(bias_calc_boxplot)

if (persist_plots) {
    bias_calc_boxplot_file <- here::here("graphs", glue("bias-boxplot-{sim_complete_time_formatted}-{sim_name}.pdf"))
    print( glue("> DEBUG: Persisting box plot to `{bias_calc_boxplot_file}`"))
    ggsave_a4(bias_calc_boxplot_file, plot = bias_calc_boxplot)
}

if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }

### Box plot of calibrated estimates using negative controls only -----

bias_calc_negctrl_boxplot <- plot_bias_calc_boxplot(
    bias_calc_sims = bias_calc_negctrl_only,
    plot_title = "Boxplot of biases (using negative controls only)",
    sim_desc = common_plot_subtitle,
    axis_break_gap = 0.025)

bias_calc_negctrl_boxplot_zoomed <- zoom_in_bias_calc_boxplot(
    bias_boxplot = bias_calc_negctrl_boxplot,
    bias_calc_sims = bias_calc_negctrl_only,
    quantile_prob_limits = c(0.0, 0.95),
    bias_colname = "est_bias",
    plot_title = "Boxplot of biases (using negative controls only) zoomed 0.0-0.95 quantiles",
    sim_desc = common_plot_subtitle)

print(bias_calc_negctrl_boxplot_zoomed)

if (persist_plots) {
    bias_calc_negctrl_boxplot_file <- here::here("graphs",
        glue("bias_negctrlonly-boxplot-{sim_complete_time_formatted}-{sim_name}.pdf"))
    print( glue("> DEBUG: Persisting box plot to `{bias_calc_negctrl_boxplot_file}`"))

    ggsave_a4(bias_calc_negctrl_boxplot_file, plot = bias_calc_negctrl_boxplot_zoomed)
}

if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }


## Histogram of biases ----
bias_calc_hist <- plot_bias_calc_histogram(bias_calc_sims = bias_calc_from_ctrls,
                                           sim_desc = common_plot_subtitle)
print(bias_calc_hist)
if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }

if (persist_plots) {
    bias_calc_hist_file <- here::here("graphs", glue("bias-histogram-{sim_complete_time_formatted}-{sim_name}.pdf"))
    print( glue("> DEBUG: Persisting histogram plot to `{bias_calc_hist_file}`"))
    ggsave_a4(bias_calc_hist_file, plot = bias_calc_hist)
}

## Density of biases - better for large simulations ----
bias_density_plot <- plot_bias_calc_density(bias_calc_sims = bias_calc_from_ctrls,
                                           sim_desc = common_plot_subtitle)
print(bias_density_plot)
if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }

if (persist_plots) {
    bias_calc_density_file <- here::here("graphs", glue("bias-density-{sim_complete_time_formatted}-{sim_name}.pdf"))
    print( glue("> DEBUG: Persisting density plot to `{bias_calc_density_file}`"))
    ggsave_a4(bias_calc_density_file, plot = bias_density_plot)
}

print("> DEBUG: Creating funnel plot")

## Funnel plot showing true estimates at different CIs ----
trt_est_from_sim <- bias_calc_from_ctrls %>% tidyr::unnest(cols = c(res))
funnel_plot <- trt_est_from_sim %>% plot_95ci_coverage_funnel(
    estimates_var = "est", se_var = "est_stderr",
    type_var = "restype",
    true_est = trt_eff,
    max_se_limit_quantile_prob = 0.95,
    plot_title = "Funnel plot of estimates and their standard errors zoomed in at 0.95 quantile of SErr",
    plot_subtitle = common_plot_subtitle)
print(funnel_plot)

if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }

if (persist_plots) {
    funnel_plot_file <- here::here("graphs", glue("funnel-plot-{sim_complete_time_formatted}-{sim_name}.pdf"))
    print( glue("> DEBUG: Persisting funnel plot to `{funnel_plot_file}`"))
    ggsave_a4(funnel_plot_file, plot = funnel_plot)
}

## Funnel plot showing BIASES and SEs of estimates ----
bias_v_se_funnelplot <- plot_bias_v_estse(combined_bias_est_data = calib_result_from_sim,
                              max_se_quantile = 0.95)
subtitle_text <- count_na_bias_calc_tbl(calib_result_from_sim)

bias_v_se_funnelplot <- bias_v_se_funnelplot +
       ggtitle("Funnel plot of bias vs SE of estimates. SE limited to 95th quantile",
        subtitle = subtitle_text)
print(bias_v_se_funnelplot)
if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }

if (persist_plots) {
    funnel_plot_file_2 <- here::here("graphs", glue("funnel-plot-bias_v_se-{sim_complete_time_formatted}-{sim_name}.pdf"))
    print( glue("> DEBUG: Persisting bias vs estimate SE funnel plot to `{funnel_plot_file_2}`"))
    ggsave_a4(funnel_plot_file_2, plot = bias_v_se_funnelplot)
}

## Empirical calibration error model plots ----

print("> DEBUG: Creating error model plot")
ctrl_est_tbl <- extract_tbl_from_sim_result(sim_res = sim_out, name = "all_ctrl_est_tbl")
ctrl_est_plot_data <- ctrl_est_tbl[['res']][[sim_num_to_plot]]

sys_err_model_plot <- EmpiricalCalibration::plotErrorModel(
    logRr = ctrl_est_plot_data[['actual_eff_est']],
    seLogRr = ctrl_est_plot_data[['actual_eff_est_stderr']],
    trueLogRr = ctrl_est_plot_data[['true_eff']]) +
    labs(subtitle = glue("Sim #{sim_num_to_plot}"))
print(sys_err_model_plot)
if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }

ctrl_est_forest_plot <- plotTrueAndObserved(ctrl_est_plot_data[['actual_eff_est']],
                                            ctrl_est_plot_data[['actual_eff_est_stderr']],
                                            ctrl_est_plot_data[['true_eff']])
ctrl_est_forest_plot <- ctrl_est_forest_plot + xlab("True effect size") +
    ggtitle("Forest plot of effect sizes from controls",
            subtitle = glue("Sim #{sim_num_to_plot}. Orange: Statistically significant estimate"))
print(ctrl_est_forest_plot)
if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }

## Funnel plot of calibrated controls ----

### Funnel plot of calibrated negative controls ----

negctrl_only_syserr_model <- extract_tbl_from_sim_result(
    sim_res = sim_out, name = "sys_err_model_negctrl_only")
null_syserr_model <- extract_tbl_from_sim_result(
    sim_res = sim_out, name = "sys_err_model_from_null")
adjtrue_posctrl_syserr_model <- extract_tbl_from_sim_result(
    sim_res = sim_out, name = "sys_err_model_adj_true_eff")
allctrl_syserr_model <- extract_tbl_from_sim_result(
    sim_res = sim_out, name = "sys_err_allctrls")

#### Funnel plot of calibrated negative controls: Using neg ctrl only sys err model ----
print("> DEBUG: Creating funnel plot of calibrated negative controls: Using neg ctrl only sys err model")
negctrl_calib_res_negctrlonly <- extract_tbl_from_sim_result(
    sim_res = sim_out, name = "calibrated_negctrl_est_negctrlonlysyserr")
negctrl_calib_funnel_plot <- plot_calib_n_uncalib_control_empclib_funnelplots(
    negctrl_calib_res_negctrlonly[['res']][[sim_num_to_plot]],
    syserr_model = negctrl_only_syserr_model[['res']][[sim_num_to_plot]],
    syserr_model_uselegacy = sys_err_model_legacy_flag,
    title = glue("Sim# {sim_num_to_plot}: NegCtrl calibrated using Neg Controls only sys err model"))
print(negctrl_calib_funnel_plot)
if (persist_plots) {
    negcalib_plotfile <- here::here("graphs", glue("funnel_plot-calib_negctrl_negctrlonly-{sim_complete_time_formatted}-{sim_name}.pdf"))
    print( glue("> DEBUG: Persisting calibrated negative control using neg ctrls only to `{negcalib_plotfile}`"))
    ggsave_a4r(negcalib_plotfile, plot = negctrl_calib_funnel_plot)
}
if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }

#### Funnel plot of calibrated negative controls: Using NULL sys err model ----
print("> DEBUG: Creating funnel plot of calibrated negative controls: Using NULL sys err model")
negctrl_calib_nullmodel <- extract_tbl_from_sim_result(
    sim_res = sim_out, name = "calibrated_negctrl_est_nullsyserr")
negctrl_calib_nullmodel_funnel_plot <- plot_calib_n_uncalib_control_empclib_funnelplots(
    negctrl_calib_nullmodel[['res']][[sim_num_to_plot]],
    syserr_model = null_syserr_model[['res']][[sim_num_to_plot]],
    syserr_model_uselegacy = sys_err_model_legacy_flag,
    title = glue("Sim# {sim_num_to_plot}: NegCtrl calibrated using Null sys err model"))
print(negctrl_calib_nullmodel_funnel_plot)
if (persist_plots) {
    negcalib_plotfile <- here::here("graphs", glue("funnel_plot-calib_negctrl_nullmodel-{sim_complete_time_formatted}-{sim_name}.pdf"))
    print( glue("> DEBUG: Persisting calibrated negative control using NULL model `{negcalib_plotfile}`"))
    ggsave_a4r(negcalib_plotfile, plot = negctrl_calib_nullmodel_funnel_plot)
}
if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }


#### Funnel plot of calibrated negative controls: Using all ctrls sys err model ----
print("> DEBUG: Creating funnel plot of calibrated negative controls: Using all ctrl sys err model")
negctrl_calib_allctrls <- extract_tbl_from_sim_result(
    sim_res = sim_out, name = "calibrated_negctrl_est")
negctrl_calib_allctrls_funnel_plot <- plot_calib_n_uncalib_control_empclib_funnelplots(
    negctrl_calib_allctrls[['res']][[sim_num_to_plot]],
    syserr_model = allctrl_syserr_model[['res']][[sim_num_to_plot]],
    syserr_model_uselegacy = sys_err_model_legacy_flag,
    title = "Sim# {sim_num_to_plot}: NegCtrl calibrated using all controls sys err model")
print(negctrl_calib_allctrls_funnel_plot)
if (persist_plots) {
    negcalib_plotfile <- here::here("graphs", glue("funnel_plot-calib_negctrl_allctrls-{sim_complete_time_formatted}-{sim_name}.pdf"))
    print( glue("> DEBUG: Persisting calibrated negative control using all controls `{negcalib_plotfile}`"))
    ggsave_a4r(negcalib_plotfile, plot = negctrl_calib_allctrls_funnel_plot)
}
if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }

### Funnel plot of calibrated positive controls ----

#### Funnel plot of positive controls: Using neg ctrl only sys err model ----
print("> DEBUG: Creating funnel plot of calibrated +positive+ controls: Using neg ctrl only sys err model")
posctrl_calib_res_negctrlonly <- extract_tbl_from_sim_result(
    sim_res = sim_out, name = "calibrated_posctrl_est_negctrlonlysyserr")
posctrl_calib_negctrlonly_funnel_plot <- plot_calib_n_uncalib_control_empclib_funnelplots(
    posctrl_calib_res_negctrlonly[['res']][[sim_num_to_plot]],
    syserr_model = negctrl_only_syserr_model[['res']][[sim_num_to_plot]],
    syserr_model_uselegacy = sys_err_model_legacy_flag,
    title = glue("Sim# {sim_num_to_plot}: PosCtrl calibrated using Neg Controls only sys err model"),
    posctrl_data = TRUE, nposctrls = length(posctrl_log_or_targets))
print(posctrl_calib_negctrlonly_funnel_plot)
if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }

if (persist_plots) {
    postctrlcalib_plotfile <- here::here("graphs", glue("funnel_plot-calib_posctrl_negctrlonly-{sim_complete_time_formatted}-{sim_name}.pdf"))
    print( glue("> DEBUG: Persisting calibrated positive control using neg ctrl only sys err model `{postctrlcalib_plotfile}`"))
    ggsave_a4r(postctrlcalib_plotfile, plot = posctrl_calib_negctrlonly_funnel_plot)
}

#### Funnel plot of positive controls: Using NULL sys err model ----
print("> DEBUG: Creating funnel plot of calibrated +positive+ controls: Using NULL sys err model")
posctrl_calib_res_nullmodel <- extract_tbl_from_sim_result(
    sim_res = sim_out, name = "calibrated_posctrl_est_nullsyserr")
posctrl_calib_nullmodel_funnel_plot <- plot_calib_n_uncalib_control_empclib_funnelplots(
    posctrl_calib_res_nullmodel[['res']][[sim_num_to_plot]],
    syserr_model = null_syserr_model[['res']][[sim_num_to_plot]],
    syserr_model_uselegacy = sys_err_model_legacy_flag,
    title = glue("Sim# {sim_num_to_plot}: PosCtrl calibrated using NULL only model"),
    posctrl_data = TRUE, nposctrls = length(posctrl_log_or_targets))
print(posctrl_calib_nullmodel_funnel_plot)
if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }

if (persist_plots) {
    postctrlcalib_plotfile <- here::here("graphs", glue("funnel_plot-calib_posctrl_nullmodel-{sim_complete_time_formatted}-{sim_name}.pdf"))
    print( glue("> DEBUG: Persisting calibrated positive control using NULL only sys err model `{postctrlcalib_plotfile}`"))
    ggsave_a4(postctrlcalib_plotfile, plot = posctrl_calib_nullmodel_funnel_plot)
}


#### Funnel plot of positive controls: Using all ctrls sys err model ----
print("> DEBUG: Creating funnel plot of calibrated +positive+ controls: Using all ctrl sys err model")
posctrl_calib_res_allctrl_model <- extract_tbl_from_sim_result(
    sim_res = sim_out, name = "calibrated_posctrl_est")
posctrl_calib_allctrl_funnel_plot <- plot_calib_n_uncalib_control_empclib_funnelplots(
    posctrl_calib_res_allctrl_model[['res']][[sim_num_to_plot]],
    syserr_model = allctrl_syserr_model[['res']][[sim_num_to_plot]],
    syserr_model_uselegacy = sys_err_model_legacy_flag,
    title = glue("Sim# {sim_num_to_plot}: PosCtrl calibrated using all controls sys err model"),
    posctrl_data = TRUE, nposctrls = length(posctrl_log_or_targets))
print(posctrl_calib_allctrl_funnel_plot)
if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }

if (persist_plots) {
    postctrlcalib_plotfile <- here::here("graphs", glue("funnel_plot-calib_posctrl_allctrls-{sim_complete_time_formatted}-{sim_name}.pdf"))
    print( glue("> DEBUG: Persisting calibrated positive control using all ctrl sys err model `{postctrlcalib_plotfile}`"))
    ggsave_a4(postctrlcalib_plotfile, plot = posctrl_calib_allctrl_funnel_plot)
}

## Funnel plot BIASES and SEs of estimates, faceted by type of calibration ----

plot_data <- calib_result_from_sim

#### Removes studies with NA standard errors. ----

subtitle_text <- count_na_bias_calc_tbl(plot_data)
no_na_calib_result_obj <- plot_data %>%
    filter_sims_without_na(na_col = "est_stderr", group_var = "sim_num")
plot_data <- no_na_calib_result_obj[['data_without_na']]

facet_labels_vec <- c("All controls", "Negative Controls only")
names(facet_labels_vec) <- c("nomod", "negctrlonlysyserr")

bias_v_se_funnelplot_faceted <- plot_bias_v_estse_v2(combined_bias_est_data = plot_data,
                                          max_se_quantile = 0.99,
                                          facet_labels_vec = facet_labels_vec)

bias_v_se_funnelplot_faceted <- bias_v_se_funnelplot_faceted %>%
    decorate_funnel_plot_with_95ci_v(plot_data = plot_data,
                                     max_se_limit_quantile_prob = 0.99)


bias_v_se_funnelplot_faceted <- bias_v_se_funnelplot_faceted +
    ggtitle("Funnel plot of bias vs SE of estimates. SE limited to 95th quantile",
            subtitle = subtitle_text)
print(bias_v_se_funnelplot_faceted)
if (manual_forward_plots) { invisible(readline(prompt="Press [enter] to continue")) }

if (persist_plots) {
    funnel_plot_file_3 <- here::here("graphs", glue("funnel-plot-bias_v_se_faceted-{sim_complete_time_formatted}-{sim_name}.pdf"))
    print( glue("> DEBUG: Persisting faceted bias vs estimate SE funnel plot to `{funnel_plot_file_2}`"))
    ggsave_a4r(funnel_plot_file_3, plot = bias_v_se_funnelplot_faceted)
}

## Funnel plot consisting of calibrated vs non-calibrated CONTROL estimates ----
if(TRUE == calibrate_negctrls) {

    if (debug_intrpt) { browser() }

    negctrl_est <- extract_tbl_from_sim_result(sim_res = sim_out, name = "negctrl_est")

    ### Data for calibration with all controls
    plot_calibrated_ctrl_wrapper(sim_res_data = sim_out, sim_res_name = "calibrated_negctrl_est",
                                 plot_title = "Calibrated and Uncalibrated Negative Controls (calibrated using all controls)",
                                 persist_data = TRUE,
                                 data_path = here::here("data", "sim-binaryoutcome"),
                                 data_filename = glue("negctrl_calib-{sim_complete_time_formatted}-{sim_name}.csv"),
                                 persist_plot = TRUE,
                                 plot_path = here::here("graphs"),
                                 plot_filename = glue("funnel_plot-negctrl_calib_allctrls-{sim_complete_time_formatted}-{sim_name}.pdf"),
                                 forest_plot_filename = glue("forest_plot-negctrl_calib_allctrls-{sim_complete_time_formatted}-{sim_name}.pdf"))

    ### Data for calibration using negative controls only
    plot_calibrated_ctrl_wrapper(sim_res_data = sim_out, sim_res_name = "calibrated_negctrl_est_negctrlonlysyserr",
                                 plot_title = "Calibrated and Uncalibrated Negative Controls (calibrated using negative controls only)",
                                 persist_data = TRUE,
                                 data_path = here::here("data", "sim-binaryoutcome"),
                                 data_filename = glue("negctrl_calib_negctrlonly-{sim_complete_time_formatted}-{sim_name}.csv"),
                                 persist_plot = TRUE,
                                 plot_path = here::here("graphs"),
                                 plot_filename = glue("funnel_plot-negctrl_calib_negctrlonly-{sim_complete_time_formatted}-{sim_name}.pdf"),
                                 forest_plot_filename = glue("forest_plot-negctrl_calib_negctrlonly-{sim_complete_time_formatted}-{sim_name}.pdf"))


    ### Data for calibration using adjusted truth in positive controls only
    plot_calibrated_ctrl_wrapper(sim_res_data = sim_out, sim_res_name = "calibrated_negctrl_est_adjtruesyserr",
                                 plot_title = "Calibrated and Uncalibrated Negative Controls (calibrated using adjusted true eff in positive controls)",
                                 persist_data = TRUE,
                                 data_path = here::here("data", "sim-binaryoutcome"),
                                 data_filename = glue("negctrl_calib_adjtruesyserr-{sim_complete_time_formatted}-{sim_name}.csv"),
                                 persist_plot = TRUE,
                                 plot_path = here::here("graphs"),
                                 plot_filename = glue("funnel_plot-negctrl_calib_adjtruesyserr-{sim_complete_time_formatted}-{sim_name}.pdf"),
                                 forest_plot_filename = glue("forest_plot-negctrl_calib_adjtruesyserr-{sim_complete_time_formatted}-{sim_name}.pdf"))

    if (debug_intrpt) { browser() }
} # End (TRUE == calibrate_negctrls)

if (TRUE == calibrate_posctrls & FALSE) {
    ## Calibrated vs non-calibrated +POSITIVE+ control estimates ----
    calibrated_posctrl_est <- extract_tbl_from_sim_result(sim_res = sim_out, name = "calibrated_posctrl_est")


    plot_calibrated_ctrl_wrapper(sim_res_data = sim_out, sim_res_name = "calibrated_posctrl_est",
                                 plot_title = "Calibrated and Uncalibrated Positive Controls (calibrated with all controls)",
                                 persist_data = TRUE,
                                 data_path = here::here("data", "sim-binaryoutcome"),
                                 data_filename = glue("posctrl_calib_allctrls-{sim_complete_time_formatted}-{sim_name}.csv"),
                                 persist_plot = TRUE,
                                 plot_path = here::here("graphs"),
                                 plot_filename = glue("funnel_plot-posctrl_calib_allctrls-{sim_complete_time_formatted}-{sim_name}.pdf"),
                                 forest_plot_filename = glue("forest_plot-posctrl_calib_allctrls-{sim_complete_time_formatted}-{sim_name}.pdf"))

    #browser()

    plot_calibrated_ctrl_wrapper(sim_res_data = sim_out, sim_res_name = "calibrated_posctrl_est_negctrlonlysyserr",
                                 plot_title = "Calibrated and Uncalibrated Positive Controls (calibrated using negative control only)",
                                 persist_data = TRUE,
                                 data_path = here::here("data", "sim-binaryoutcome"),
                                 data_filename = glue("posctrl_calib_negctrlonly-{sim_complete_time_formatted}-{sim_name}.csv"),
                                 persist_plot = TRUE,
                                 plot_path = here::here("graphs"),
                                 plot_filename = glue("funnel_plot-posctrl_calib_negctrlonly-{sim_complete_time_formatted}-{sim_name}.pdf"),
                                 forest_plot_filename = glue("forest_plot-posctrl_calib_negctrlonly-{sim_complete_time_formatted}-{sim_name}.pdf"))

    plot_calibrated_ctrl_wrapper(sim_res_data = sim_out, sim_res_name = "calibrated_posctrl_est_adjtruesyserr",
                                 plot_title = "Calibrated and Uncalibrated Positive Controls (calibrated using adjusted true effect)",
                                 persist_data = TRUE,
                                 data_path = here::here("data", "sim-binaryoutcome"),
                                 data_filename = glue("posctrl_calib_adjtruesyserr-{sim_complete_time_formatted}-{sim_name}.csv"),
                                 persist_plot = TRUE,
                                 plot_path = here::here("graphs"),
                                 plot_filename = glue("funnel_plot-posctrl_calib_adjtruesyserr-{sim_complete_time_formatted}-{sim_name}.pdf"),
                                 forest_plot_filename = glue("forest_plot-posctrl_calib_adjtruesyserr-{sim_complete_time_formatted}-{sim_name}.pdf"))


    ## Merging
    print("> DEBUG: Merging calibrated positive control data")
    calibrated_posctrl_allctrls <- extract_tbl_from_sim_result(sim_res = sim_out, name = "calibrated_posctrl_est") %>% unnest(c(sim_num, res))
    calibrated_posctrl_negctrlonly <- extract_tbl_from_sim_result(sim_res = sim_out, name = "calibrated_posctrl_est_negctrlonlysyserr") %>% unnest(c(sim_num, res))
    calibrated_posctrl_adjtruesyserr <- extract_tbl_from_sim_result(sim_res = sim_out, name = "calibrated_posctrl_est_adjtruesyserr") %>% unnest(c(sim_num, res))
    posctrl_calib_merged_data <- rbind(calibrated_posctrl_allctrls , calibrated_posctrl_negctrlonly  , calibrated_posctrl_adjtruesyserr)
    readr::write_csv(posctrl_calib_merged_data,
                     here::here("data", "sim-binaryoutcome", glue("posctrl_calibrated-{sim_complete_time_formatted}-{sim_name}.csv")))

} # End (TRUE == calibrate_posctrls)
