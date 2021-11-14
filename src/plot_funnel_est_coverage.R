#' Funnel plot of estimates and their 95% confidence interval
#'
#' See
#' https://sakaluk.wordpress.com/2016/02/16/7-make-it-pretty-plots-for-meta-analysis/
#'https://stats.stackexchange.com/a/195333
#' @export
plot_95ci_coverage_funnel <- function(data, estimates_var = "estimate", se_var = "se",
                                      type_var = "type",
                                      true_est = NA_real_,
                                      max_se_limit_quantile_prob = 1.0,
                                      plot_title = "Funnel plot of estimates and their standard errors zoomed in at 0.95 quantile of SErr",
                                      plot_subtitle = "",
                                      num_se_funnel_points = NA_integer_) {

    if (!require(dplyr))   { library(dplyr) }
    if (!require(ggplot2)) { library(ggplot2) }
    if (!require(tibble))  { library(tibble) }

    # Sequence of standard errors to calculate the values for the 'funnel' lines

    #max_se <- data %>% filter(!is.na(!!sym(se_var))) %>% pull(!!sym(se_var)) %>%
    #    quantile(probs = c(max_se_limit_quantile_prob)) %>% .[[1]]

    max_se <- data %>% pull(!!sym(se_var)) %>% quantile(probs = c(max_se_limit_quantile_prob)) %>% .[[1]]

    if (!is.na(num_se_funnel_points)) {
        se_seq <- seq(from = 0.0, to = max_se, length.out = num_se_funnel_points)
    } else {
        se_seq <- seq(from = 0.0, to = max_se, by = 0.001)
    }
    ll95 <- true_est - (1.96 * se_seq)
    ul95 <- true_est + (1.96 * se_seq)

    ci_paths <- tibble(x = c(true_est - qnorm(0.995) * max_se,
                             true_est, true_est + qnorm(0.995) * max_se),
                       y = c(max_se, 0 , max_se))

    funnel_lines_data <- tibble(se_seq = se_seq, ll95 = ll95, ul95 = ul95)
    plot_data <- data %>% filter(!is.na(!!sym(se_var))) %>%
        filter(.data[[se_var]] <= max_se)

    funnel_plot <- plot_data %>% ggplot(aes(x = !!sym(estimates_var), y = !!sym(se_var))) +
        # geom_path(aes(x = x, y = y), data = ci_paths, linetype = "dashed") +
        geom_path(aes(x = ll95, y = se_seq), data = funnel_lines_data, linetype = "dashed") +
        geom_path(aes(x = ul95, y = se_seq), data = funnel_lines_data, linetype = "dashed") +
        geom_point(aes(colour = !!sym(type_var)), alpha = 0.7) +
        geom_vline(xintercept = true_est, linetype = "dotted") +
        scale_colour_discrete(
            name = "Type",
            #breaks = c("outcome_of_interest−calibrated", "outcome_of_interest−uncalibrated"),
            labels = c("Calibrated", "Uncalibrated")) +
        #geom_line(aes(x = se_seq, y = ll95), linetype = 'dotted', data = funnel_lines_data) +
        #geom_line(aes(x = se_seq, y = ul95), linetype = 'dotted', data = funnel_lines_data) +
        xlab("Estimates log odds ratio") + ylab("Standard Error") +
        ggtitle(label = plot_title, subtitle = plot_subtitle)

    return(funnel_plot)
}

#' Returns text on number of NAs in bias calculation.
#'
#' @param combined_bias_est_data Same data supplied to function such as
#'     \code{plot_bias_v_estse}.
count_na_bias_calc_tbl <- function(combined_bias_est_data = NULL) {
    if (!require(dplyr)) { library(dplyr) }
    if (!require(glue))  { library(glue) }

    na_count_tbl <- combined_bias_est_data %>% group_by(calitype) %>%
        summarise(nacount = sum(is.na(est_stderr)))

    return(glue(
        "No modification: {na_count_tbl %>% filter(calitype == 'nomod') %>% pull(nacount)}, ",
        "Calibrated with Negative controls only: {na_count_tbl %>% filter(calitype == 'negctrlonlysyserr') %>% pull(nacount)}, ",
        "Calibrated with adjusted true value in positive control: {na_count_tbl %>% filter(calitype == 'adjtruthposctrlsyserr') %>% pull(nacount)}"))
}

#' Plots \emph{bias} vs standard error of estimates, across types
#' of calibrations.
#'
#' @param combined_bias_est_data \code{data.frame} like data structure that
#'      contains the following columns: \emph{est_bias}, \emph{est_stderr},
#'      \emph{restype} (calibrated or uncalibrated estimate, mapped to colours),
#'      \emph{calitype} (calibrated with all controls, negative controls only,
#'      or with adjusted true values of positive controls, mapped to
#'      shape of points).
#'
plot_bias_v_estse <- function(combined_bias_est_data = NULL,
                              max_se_quantile = 0.95) {

    if (!require(ggplot2)) { library(ggplot2) }
    if (!require(glue))    { library(glue) }

    se_95quantile_na_removed <- quantile(combined_bias_est_data[['est_stderr']] ,
                                         max_se_quantile, na.rm = TRUE)

    ggplot(data = combined_bias_est_data, aes(x = est_bias, y = est_stderr)) +
        geom_point(aes(colour = restype , shape = calitype), alpha = 0.8, size = 1.7) +
        geom_vline(xintercept = 0.0) +
        scale_shape_manual(name = "Calibration Type",
                           labels = c("AllControls", "NegativeControlsOnly", "AdjustedTruthPositiveControl"),
                           breaks = c("nomod", "negctrlonlysyserr", "adjtruthposctrlsyserr"),
                           values = c(2, 3, 4)) +
        scale_color_discrete(name = "Estimate type") +
        coord_cartesian(ylim = c(0, se_95quantile_na_removed)) +
        xlab("Bias") +
        ylab("Standard Error of estimate")
}

#' Plots \emph{bias} vs standard error of estimates, \emph{faceted} by types
#' of calibrations.
#'
#' @param combined_bias_est_data \code{data.frame} like data structure that
#'      contains the following columns: \emph{est_bias}, \emph{est_stderr},
#'      \emph{restype} (calibrated or uncalibrated estimate, mapped to colours),
#'      \emph{calitype} (calibrated with all controls, negative controls only,
#'      or with adjusted true values of positive controls, mapped to
#'      shape of points).
#' @param max_se_quantile Maximum quantile of standard error to plot. Defaults to
#'      95%, represented as 0.95
#' @param facet_labels_vec Named vector specifying the labels of the facet.
#'      The 'names' component of the vector corresponds to the values of the facet variable.
#'      The values component contains the values to display.
plot_bias_v_estse_v2 <- function(combined_bias_est_data = NULL,
                                 max_se_quantile = 0.95,
                                 facet_labels_vec = c(),
                                 zoom_max_se_quantils = TRUE) {

    if (!require(ggplot2)) { library(ggplot2) }
    if (!require(glue))    { library(glue) }

    se_95quantile_na_removed <- quantile(combined_bias_est_data[['est_stderr']] ,
                                         max_se_quantile, na.rm = TRUE)

    bias_v_estse_plot <- ggplot(data = combined_bias_est_data, aes(x = est_bias, y = est_stderr)) +
        geom_point(aes(colour = restype), alpha = 0.8, size = 1.7) +
        geom_vline(xintercept = 0.0) +
        xlab("Bias") +
        ylab("Standard Error of estimate")
        
    if (zoom_max_se_quantils) { 
        bias_v_estse_plot <- bias_v_estse_plot +
		 coord_cartesian(ylim = c(0, se_95quantile_na_removed))
    }

    if (length(facet_labels_vec) > 0 & !is.null(facet_labels_vec)) {
        bias_v_estse_plot <- bias_v_estse_plot + facet_wrap(vars(calitype),
                                                            labeller = labeller(calitype = facet_labels_vec) )
    } else {
        bias_v_estse_plot <- bias_v_estse_plot + facet_wrap(vars(calitype))
    }
    return(bias_v_estse_plot)

} # End plot_bias_v_estse_v2()

#' Decorates bias vs estimate standard error plot with 95 CI lines.
#'
#' @param funnelplot Bias vs Estimate Standard Error funnel plot
#' @param plot_data Data used to generated \code{funnelplot}.
#' @param max_se_limit_quantile_prob Max quantile cut off of standard error
#'
#' @return \code{ggplot} object with 95CI V-line
decorate_funnel_plot_with_95ci_v <- function(funnelplot = NULL,
                                             plot_data = NULL, manual_max_se = NULL,
                                             max_se_limit_quantile_prob = 0.99) {

    if (is.null(manual_max_se)) {
        if (is.null(plot_data)) { return(NULL) }
        # Calculates required v-line data , taking into account specified quantile ----
        max_se <- plot_data %>% pull(est_stderr) %>% quantile(probs = c(max_se_limit_quantile_prob)) %>% .[[1]]
    } else {
        # Uses user supplied level of maximum value of standard error, ignores plot data and quantiles.
        max_se <- manual_max_se
    }

    bias_centre_value = 0.0
    se_seq <- seq(from = 0.0, to = max_se, by = 0.01)
    ll95 <- bias_centre_value - (1.96 * se_seq)
    ul95 <- bias_centre_value + (1.96 * se_seq)

    #ci_paths <- tibble(x = c(bias_centre_value - qnorm(0.95) * max_se,
    #                         bias_centre_value,
    #                         bias_centre_value + qnorm(0.95) * max_se),
    #                   y = c(max_se, 0 , max_se))

    funnel_lines_data <- tibble(se_seq = se_seq, ll95 = ll95, ul95 = ul95)

    # Applies v-line decoration ----
    plot_to_rt <- funnelplot +
        geom_path(aes(x = ll95, y = se_seq), data = funnel_lines_data, linetype = "dashed") +
        geom_path(aes(x = ul95, y = se_seq), data = funnel_lines_data, linetype = "dashed")
} # End decorate_funnel_plot_with_95ci_v()

#' Plots funnel plot of calibrated control
#'
plot_calibrated_ctrl_wrapper <- function(sim_res_data = NULL, sim_res_name = "",
                                         plot_title = "",
                                         persist_data = FALSE,
                                         data_path = here::here("data"),
                                         data_filename = "data.csv",
                                         persist_plot = FALSE,
                                         plot_path = here::here("graphs"),
                                         plot_filename = "funnelplot.pdf",
                                         forest_plot_filename = "forestplot.pdf") {

    if (!require(readr)) { library(readr) }

    calibrated_negctrl_est <- extract_tbl_from_sim_result(sim_res = sim_res_data, name = sim_res_name)

    negctrl_calib_plotting_tbl_raw <- calibrated_negctrl_est %>% tidyr::unnest(c(sim_num, res))
    if (TRUE == persist_data) {
        negctrl_calib_data_file <- file.path(data_path, data_filename)
        print( glue("> DEBUG: Persiting calibrated -negctrl- data to `{negctrl_calib_data_file}`"))
        readr::write_csv(negctrl_calib_plotting_tbl_raw, negctrl_calib_data_file)
    }

    # Removes studies with NA standard errors
    no_na_data_obj <- negctrl_calib_plotting_tbl_raw %>%
        filter_sims_without_na(na_col = "actual_eff_est_stderr", group_var = "sim_num")
    num_na <- no_na_data_obj[['na_groups']] %>% length()

    plot_subtitle <- glue("Across all negative control outcomes, and all simulations\n",
                          "{num_na} simulations with NA SE removed")

    ### Funnel plot calibrated negative controls ----
    print("> DEBUG: Creating funnel plot of calibrated negative controls")

    negctrl_calib_plotting_tbl <- no_na_data_obj[['data_without_na']]
    negctrl_calib_res_funnel_plot <- negctrl_calib_plotting_tbl %>% plot_95ci_coverage_funnel(
        estimates_var = "actual_eff_est", se_var = "actual_eff_est_stderr",
        type_var = "restype",
        true_est = 0.0,
        num_se_funnel_points = nrow(negctrl_calib_plotting_tbl))

    negctrl_calib_res_funnel_plot <- negctrl_calib_res_funnel_plot +
        ggtitle(plot_title, subtitle = plot_subtitle)
    print(negctrl_calib_res_funnel_plot)

    if (persist_plot) {
        negctrl_funnel_plot_file <- file.path(plot_path, plot_filename)
        print( glue("> DEBUG: Persisting negative control calibration funnel plot to `{negctrl_funnel_plot_file}`"))
        ggsave_a4(negctrl_funnel_plot_file, plot = negctrl_calib_res_funnel_plot)
    }

    ### Forest plot of calibrated negative controls ----

    #negctrl_calib_plot_data_without_na_se <- negctrl_calib_plotting_tbl %>%
    #    filter_sims_without_na(na_col = "actual_eff_est_stderr", group_var = "sim_num")
    #negctrl_calib_plot_data_na_se_count <- negctrl_calib_plot_data_without_na_se[['na_groups']] %>%
    #    length()

    negctrl_calib_forest_plot <- negctrl_calib_plotting_tbl %>% plot_forest_general(
        estimate_colname = "actual_eff_est", stderr_colname = "actual_eff_est_stderr",
        lower95_ci_colname = "est_L95", upper95_ci_colname = "est_U95",
        aes_arrange_var = "restype",
        true_est_val = 0.0,
        plot_title = plot_title,
        plot_subtitle = plot_subtitle)
    print(negctrl_calib_forest_plot)

    if (persist_plot) {
        negctrl_forest_plot_file <- file.path(plot_path, forest_plot_filename)
        print( glue("> DEBUG: Persisting negative control calibration forest plot to `{negctrl_forest_plot_file}`"))
        ggsave_a4(negctrl_forest_plot_file, plot = negctrl_calib_forest_plot)
    }

} # End plot_calibrated_ctrl_wrapper()
