## ----------------------------------------------------------------------------
##
## Creates a CSV file that contains paths to simulation results, for
## a specified simulation scenario.
##
## ----------------------------------------------------------------------------

library(dplyr)
library(tibble)
library(glue)
library(here)
library(lubridate)
library(readr)

# Top level result directory.
data_path_toplvl <- here("data", "simresdata")

# Each sub directory under `data_path_toplvl` is named after a simulation scenario
sim_scenario <- "unmeasured_conf"   # "quadratic", "x1x2intr", "measurement_error", "nonpositivity"

negctrl_causal_types <- c("negctrl_causalstruct", "negctrl_nocausalstructure")

negctrl_nums <- c("negctrl_05", "negctrl_30")

coef_assign_cases <- c("allrand", "ideal")

if (sim_scenario != "measurement_error") {
    data_paths_df <- as.data.frame(
        expand.grid(negctrl_causal_types, negctrl_nums, coef_assign_cases))
    colnames(data_paths_df) <- c("negctrl_causal_type", "negctrl_num", "coef_assign_case")
} else if (sim_scenario == "measurement_error") {
    merror_types <- c("norand_err_gaussian_params", "rand_err_gaussian_params")

    data_paths_df <- as.data.frame(
        expand.grid(negctrl_causal_types, merror_types, negctrl_nums, coef_assign_cases))
    colnames(data_paths_df) <- c("negctrl_causal_type", "merror_type", "negctrl_num", "coef_assign_case")
}


data_paths <- apply(data_paths_df, MARGIN = 1, FUN = function(path_component, toplvlpath, scenario) {

    if (scenario == "measurement_error") {
        sim_data_path <- file.path(toplvlpath, scenario,
                                   path_component[[1]], path_component[[2]], path_component[[3]],
                                   path_component[[4]],
                                   "data")
    } else {
        sim_data_path <- file.path(toplvlpath, scenario,
                                   path_component[[1]], path_component[[2]], path_component[[3]],
                                   "data")
    }


    data_files <- list.files(path = sim_data_path, pattern = ".*bias_combined_all.*\\.csv$")
    if (length(data_files) == 0) {
        print(glue("!! WARNING: Data file not found in `{sim_data_path}`"))
        if (scenario != "measurement_error") {
            return(c(sim_scenario = scenario,
                     negctrl_causal_type = path_component[[1]],
                     negctrl_num = path_component[[2]],
                     coef_assign_case = path_component[[3]],
                     data_file_path = NA_character_))
        } else if (scenario == "measurement_error") {
            return(c(sim_scenario = scenario,
                     negctrl_causal_type = path_component[[1]],
                     merror_type = path_component[[2]],
                     negctrl_num = path_component[[3]],
                     coef_assign_case = path_component[[4]],
                     data_file_path = NA_character_))
        }
    }
    data_file <- data_files[[1]]
    data_file_path <- file.path(sim_data_path, data_file)

    path_data_to_rt <- c()
    if (scenario != "measurement_error") {
        path_data_to_rt <- c(sim_scenario = scenario,
                             negctrl_causal_type = path_component[[1]],
                             negctrl_num = path_component[[2]],
                             coef_assign_case = path_component[[3]],
                             data_file_path = data_file_path)

    } else if (scenario == "measurement_error") {
        path_data_to_rt <- c(sim_scenario = scenario,
                             negctrl_causal_type = path_component[[1]],
                             merror_type = path_component[[2]],
                             negctrl_num = path_component[[3]],
                             coef_assign_case = path_component[[4]],
                             data_file_path = data_file_path)
    }
    # browser()
    return(path_data_to_rt)

}, toplvlpath = data_path_toplvl, scenario = sim_scenario)

data_paths_fmt <- t(data_paths) %>% as_tibble()

print("> DEBUG: Path to files created")
print(head(data_paths_fmt))

print("> DEBUG: Persisting path to files")
# Location to save collated paths.
persist_dir <- here("data", "resdata-paths", lubridate::today())
if (FALSE == dir.exists(persist_dir)) {
    print(glue("> DEBUG: Creating `{persist_dir}`"))
    dir.create(persist_dir)
}

res_files_fullpath <- file.path(persist_dir, glue("res_data-paths-{sim_scenario}.csv"))
print(glue("> DEBUG: Persisting result file paths to `{res_files_fullpath}`"))
readr::write_csv(data_paths_fmt, res_files_fullpath)
