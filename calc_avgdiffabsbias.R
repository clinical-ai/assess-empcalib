# -------------------------------------------------------------------------
# Calculates average of mean standardised biases
# -------------------------------------------------------------------------

library(dplyr)
library(readr)
library(here)
library(glue)
library(purrr)
library(lubridate)

# Function definitions ----------------------------------------------------

calc_avgabsbias <- function(sim_scenario = NULL,
 	                    negctrl_causal_type = NULL,
 	                    merror_type = NA_character_,
 	                    negctrl_num = NULL,
 	                    coef_assign_case = NULL,
 	                    data_file_path = NULL) {

  print(glue("Processing {sim_scenario}-{negctrl_causal_type}-",
             "{merror_type}-{negctrl_num}-{coef_assign_case}"))

  if (is.na(data_file_path)) {
  	print("!! WARNING: No data file path, returning NA")

  	return(list(sim_scenario = sim_scenario,
  	            negctrl_causal_type = negctrl_causal_type,
  	            merror_type = merror_type,
  	            negctrl_num = negctrl_num,
                    calib_type = NA_character_,
  	            coef_assign_case = coef_assign_case,
  	            avg_std_absbias_calibrated = NA_real_,
  	            avg_std_absbias_uncalibrated = NA_real_,
  	            diff_avg_std_absbias = NA_real_))
  }

  print(glue("> Reading `{data_file_path}`"))
  simres <- read_csv(data_file_path)

  truth_colname <- "true_eff"
  est_colname <- "est"

  simres <- simres %>% mutate(std_abs_bias =
    abs(simres[[truth_colname]] - simres[[est_colname]]) /
    simres[[truth_colname]])

  calc_diff_avg_bias <- function(simres = NULL) {
    print("> DEBUG: calc_diff_avg_bias() called")

    res_calibration_colname <- "restype"
    calibrated_res <- simres %>%
      filter(.data[[res_calibration_colname]] == "calibrated")
    uncalibrated_res <- simres %>%
      filter(.data[[res_calibration_colname]] == "uncalibrated")

    avg_std_abs_bias_calibrated <- simres %>%
      filter(.data[[res_calibration_colname]] == "calibrated") %>%
      pull(std_abs_bias) %>% mean()

    avg_std_abs_bias_uncalibrated <- simres %>%
      filter(.data[[res_calibration_colname]] == "uncalibrated") %>%
      pull(std_abs_bias) %>% mean()

    diff_avg_std_abs_bias <- avg_std_abs_bias_calibrated - avg_std_abs_bias_uncalibrated

    return(list(avg_bias_calibrated = avg_std_abs_bias_calibrated,
                avg_bias_uncalibrated = avg_std_abs_bias_uncalibrated,
                diff_avg_abs_bias = diff_avg_std_abs_bias))
  }

  syserrmodel_colname <- "calitype"

  null_syserrmodel_res <- simres %>%
    filter(.data[[syserrmodel_colname]] == "negctrlonlysyserr") %>%
    calc_diff_avg_bias()

  allctrl_syserrmodel_res <- simres %>%
    filter(.data[[syserrmodel_colname]] == "nomod") %>%
    calc_diff_avg_bias()

  null_syserrmodel_restbl <- tibble(
    sim_scenario = sim_scenario,
    negctrl_causal_type = negctrl_causal_type,
    merror_type = merror_type,
    negctrl_num = negctrl_num,
    syserrmodel = "negctrlonlysyserr",
    coef_assign_case = coef_assign_case,
    avg_std_absbias_calibrated = null_syserrmodel_res[["avg_bias_calibrated"]],
    avg_std_absbias_uncalibrated = null_syserrmodel_res[["avg_bias_uncalibrated"]],
    diff_avg_std_absbias = null_syserrmodel_res[["diff_avg_abs_bias"]]
  )
  print("> Results for negative control only calibration calculated")

  allctrl_syserrmodel_restbl <- tibble(
    sim_scenario = sim_scenario,
    negctrl_causal_type = negctrl_causal_type,
    merror_type = merror_type,
    negctrl_num = negctrl_num,
    syserrmodel = "nomod",
    coef_assign_case = coef_assign_case,
    avg_std_absbias_calibrated = allctrl_syserrmodel_res[["avg_bias_calibrated"]],
    avg_std_absbias_uncalibrated = allctrl_syserrmodel_res[["avg_bias_uncalibrated"]],
    diff_avg_std_absbias = allctrl_syserrmodel_res[["diff_avg_abs_bias"]]
  )
  print("> Results for all controls calibration calculated")

  res_tbl <- rbind(null_syserrmodel_restbl, allctrl_syserrmodel_restbl)
  return(res_tbl)
}

# Main --------------------------------------------------------------------

simres_filepath_info_file <- here("data", "resdata-paths",
  "res_data-paths-all_scenarios.csv")

filepath_info_file <- simres_filepath_info_file
print(glue("Reading paths from `{filepath_info_file}`"))
file_paths <- read_csv(filepath_info_file)

scenario_to_process <- NA_character_
if (!is.na(scenario_to_process)) {
  file_paths <- file_paths %>% filter(sim_scenario == scenario_to_process)
}
file_paths <- file_paths %>% arrange(negctrl_causal_type, negctrl_num, coef_assign_case)

calc_result <- file_paths %>% pmap_dfr(calc_avgabsbias)
print(calc_result)

result_file <- "simres_avg_diffabsbias_tbl2_data.csv"
result_persist_dir <- here("data", "simres_avg_diffabsbias", lubridate::today())
if (!dir.exists(result_persist_dir)) {
    print(glue("Creating directory `{result_persist_dir}`"))
    dir.create(result_persist_dir, recursive = TRUE)
}
persist_path <- file.path(result_persist_dir, result_file)
print(glue("> Persisting result to `{persist_path}`"))
write_csv(calc_result, persist_path)

