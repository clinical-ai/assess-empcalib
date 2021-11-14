## ----------------------------------------------------------------------------
##
## Code to generate funnel plots figures with a custom value for
## y-axis (standard error) value.
##
## ----------------------------------------------------------------------------

library(dplyr)
library(tibble)
library(readr)
library(glue)
library(here)

library(ggplot2)

source(here("src", "plot_funnel_est_coverage.R"))

# Find locations of collated result files -----
alldata_respath_file <- here("data/resdata-paths/res_data-paths-all_scenarios.csv")
print(glue("> Reading result files from `{alldata_respath_file}`"))
simres_filepaths <- read_csv(alldata_respath_file)

## Specify persist location -----
graph_loc <- here("graphs", "bias_v_se-manual_yaxis")

# Declare functions -----
gen_new_funnel_plot <- function(plotdata = NULL, manual_max_se = 0.0, drop_na_stderr = FALSE) {
    if (drop_na_stderr) {
        plotdata <- plotdata %>% tidyr::drop_na(est_stderr)
    }
    print("> Generating initial funnel plot")
    funnelplot <- plot_bias_v_estse_v2(
        plotdata,
        max_se_quantile = 0.95,
        zoom_max_se_quantils = FALSE,
        facet_labels_vec = c('nomod' = 'All controls', 'negctrlonlysyserr' = 'Negative Controls only')) +
        ylim(c(0.0, manual_max_se)) +
        ggtitle(label = "Funnel plot of bias vs SE of estimates. SE limited to 95th quantile")

    print("> Decorating funnel plot with V-shape")
    funnelplot_with_v <- funnelplot %>%
        decorate_funnel_plot_with_95ci_v(manual_max_se = manual_max_se)

    return(funnelplot_with_v)
} # End gen_new_funnel_plot()

# Re-generate plots -----
## Unmeasured confounder scenario -----
## Change the filter parameters to plot the required results.

unmeasuredconf_res_file <- simres_filepaths %>%
    filter(sim_scenario == "unmeasured_conf") %>%
    filter(negctrl_causal_type == "negctrl_nocausalstructure" &  # "negctrl_causalstruct" &
           negctrl_num         == "negctrl_05" &                 # "negctrl_30" &
           coef_assign_case    == "allrand") %>%                 # "ideal") %>%
    pull(data_file_path)
unmeasuredconf_unsuitable_negctrl_5_allrand_data <- read_csv(unmeasuredconf_res_file)

print("> Generating new funnel plot for Unmeasure confounding scenario")
unmeasuredconf_unsuitable_negctrl_5_allrand_funnelplot_with_v <-
    gen_new_funnel_plot(unmeasuredconf_unsuitable_negctrl_5_allrand_data,
                        manual_max_se = 0.4,    # Limit of y-axis
                        drop_na_stderr = TRUE)

print(unmeasuredconf_unsuitable_negctrl_5_allrand_funnelplot_with_v)


unmeasured_conf_plot_filename <- "funnel_plot-bias_se-unmeasured_conf-unsuitable_negctrl-5_negctrls-allrand"

unmeasured_conf_plot_pdf_filename <- paste0(unmeasured_conf_plot_filename, ".pdf")
unmeasured_conf_plot_pdf_path <- file.path(graph_loc, unmeasured_conf_plot_pdf_filename)
print(glue("Persisting unmeasured confounder funnel plot as PDF to `{unmeasured_conf_plot_pdf_path}`"))
ggsave(unmeasured_conf_plot_pdf_path,
       plot = unmeasuredconf_unsuitable_negctrl_5_allrand_funnelplot_with_v,
       width = 297, height = 210, units = "mm") # A4 landscape

## Model-mis-spec quadratic term scenario -----
## Repeating the same process as unmeasured confounders 

quadratic_sim_res_file <- simres_filepaths %>%
    filter(sim_scenario == "quadratic") %>%
    filter(negctrl_causal_type  == "negctrl_nocausalstructure" & # "negctrl_causalstruct" &
               negctrl_num      == "negctrl_05" &                # "negctrl_30" &
               coef_assign_case == "allrand") %>%                # "ideal") %>%
    pull(data_file_path)
quadratic_sim_unsuitable_negctrl_5_allrand_data <- read_csv(quadratic_sim_res_file)

print("> Generating new funnel plot for quadratic term scenario")
quadratic_sim_unsuitable_negctrl_5_allrand_funnelplot_with_v <-
    gen_new_funnel_plot(quadratic_sim_unsuitable_negctrl_5_allrand_data,
                        manual_max_se = 0.4) # Limit of y-axis

print(quadratic_sim_unsuitable_negctrl_5_allrand_funnelplot_with_v)

quadratic_sim_plot_filename <- "funnel_plot-bias_se-quadratic-unsuitable_negctrl-5_negctrls-allrand"

quadratic_sim_plot_pdf_filename <- paste0(quadratic_sim_plot_filename, ".pdf")
quadratic_sim_plot_pdf_path <- file.path(graph_loc, quadratic_sim_plot_pdf_filename)
print(glue("Persisting quadratic term funnel plot as PDF to `{quadratic_sim_plot_pdf_path}`"))
ggsave(filename = quadratic_sim_plot_pdf_path,
       plot = quadratic_sim_unsuitable_negctrl_5_allrand_funnelplot_with_v,
       width = 297, height = 210, units = "mm")

## Model-mis-spec interaction term scenario -----

x1x2_sim_res_file <- simres_filepaths %>%
    filter(sim_scenario == "x1x2intr") %>%
    filter(negctrl_causal_type  == "negctrl_nocausalstructure" & # "negctrl_causalstruct" &
               negctrl_num      == "negctrl_05" &                # "negctrl_30" &
               coef_assign_case == "allrand") %>%                # "ideal") %>%
    pull(data_file_path)
x1x2intr_sim_unsuitable_negctrl_5_allrand_data <- read_csv(x1x2_sim_res_file)

print("> Generating new funnel plot for interaction terms scenario")
x1x2intr_sim_unsuitable_negctrl_5_allrand_funnelplot_with_v <-
    gen_new_funnel_plot(x1x2intr_sim_unsuitable_negctrl_5_allrand_data,
                        manual_max_se = 0.25,  # Limit of y-axis
                        drop_na_stderr = TRUE)
print(x1x2intr_sim_unsuitable_negctrl_5_allrand_funnelplot_with_v)

x1x2intr_sim_plot_filename <- "funnel_plot-bias_se-x1x2intr-unsuitable_negctrl-5_negctrls-allrand"

x1x2intr_sim_plot_pdf_filename <- paste0(x1x2intr_sim_plot_filename, ".pdf")
x1x2intr_sim_plot_pdf_path <- file.path(graph_loc, x1x2intr_sim_plot_pdf_filename)
print(glue("Persisting x1x2 interaction funnel plot as PDF to `{x1x2intr_sim_plot_pdf_path}`"))
ggsave(filename = x1x2intr_sim_plot_pdf_path,
       plot = x1x2intr_sim_unsuitable_negctrl_5_allrand_funnelplot_with_v,
       width = 297, height = 210, units = "mm")

## Lack of positivity term scenario -----

nonpositivity_sim_res_file <- simres_filepaths %>%
    filter(sim_scenario == "nonpositivity") %>%
    filter(negctrl_causal_type  == "negctrl_nocausalstructure" & # "negctrl_causalstruct" &
               negctrl_num      == "negctrl_05" &                # "negctrl_30" &
               coef_assign_case == "allrand") %>%                # "ideal") %>%
    pull(data_file_path)
nonpositivity_sim_unsuitable_negctrl_5_allrand_data <- read_csv(nonpositivity_sim_res_file)

print("> Generating new funnel plot for non-positivity scenario")
nonpositivity_sim_unsuitable_negctrl_5_allrand_funnelplot_with_v <-
    gen_new_funnel_plot(nonpositivity_sim_unsuitable_negctrl_5_allrand_data,
                        manual_max_se = 0.08,  # Limit of y-axis
                        drop_na_stderr = TRUE)
print(nonpositivity_sim_unsuitable_negctrl_5_allrand_funnelplot_with_v)

nonpositivity_sim_plot_filename <- "funnel_plot-bias_se-nonpositivity-unsuitable_negctrl-5_negctrls-allrand"
nonpositivity_sim_plot_filename <- paste0(nonpositivity_sim_plot_filename, ".pdf")
nonpositivity_sim_plot_pdf_path <- file.path(graph_loc, nonpositivity_sim_plot_filename)
print(glue("Persisting non-positivy interaction funnel plot as PDF to `{nonpositivity_sim_plot_pdf_path}`"))
ggsave(filename = nonpositivity_sim_plot_pdf_path,
       plot = nonpositivity_sim_unsuitable_negctrl_5_allrand_funnelplot_with_v,
       width = 297, height = 210, units = "mm")


## Measurement error term scenario ----

merr_sim_res_file <- simres_filepaths %>%
    filter(sim_scenario == "measurement_error" &
           negctrl_causal_type == "negctrl_causalstruct" &
           negctrl_num         == "negctrl_05" &
           coef_assign_case    == "allrand" &
           merror_type == "rand_err_gaussian_params") %>%
    pull(data_file_path)
merr_sim_negctrl_5_allrand_data <- read_csv(merr_sim_res_file)

print("> Generating new funnel plot for measurement error scenario")
merr_sim_negctrl_5_allrand_funnelplot_with_v <-
    gen_new_funnel_plot(merr_sim_negctrl_5_allrand_data,
                        manual_max_se = 0.08,     # Limit of y-axis
                        drop_na_stderr = FALSE)
print(merr_sim_negctrl_5_allrand_funnelplot_with_v)

merr_sim_plot_filename <- "merror-bias_se-nonpositivity-suitable_negctrl-5_negctrls-allrand"
merr_sim_plot_filename <- paste0(merr_sim_plot_filename, ".pdf")
merror_sim_plot_pdf_path <- file.path(graph_loc, merr_sim_plot_filename)
print(glue("Persisting measurement error plot as PDF to `{merror_sim_plot_pdf_path}`"))
ggsave(filename = merror_sim_plot_pdf_path,
       plot = merr_sim_negctrl_5_allrand_funnelplot_with_v,
       manual_max_se = 0.08,                      # Limit of y-axis
       width = 297, height = 210, units = "mm")