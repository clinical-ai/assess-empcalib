## ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
##
## Forest plot functions
##
## ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#' Plots \emph{forest plot} of estiamtes of calibrated and uncalibrated results.
#'
#' @param bias_calc_sims Table of biases obtained from calling
#'     \code{extract_tbl_from_sim_result()} function.
#'     Columns of sim_num, restype, est, est_stderr, true_eff, est_bias
#'
plot_est_forest_plot <- function(bias_calc_sims = NULL,
                                 plot_title = "Forest plot of estimate of treatment effect",
                                 sim_desc = "sim-desc",
                                 remove_sim_with_na_calib_stderr = FALSE,
                                 calibrated_max_est_se_quantile = 1.0) {

    trt_est_from_sim <- bias_calc_sims %>% tidyr::unnest(cols = c(res))

    num_calib_est_with_na_stderr <- NA_integer_

    trt_est_to_plot <- trt_est_from_sim
    sim_subtile <- sim_desc

    if (remove_sim_with_na_calib_stderr) {
        na_stderr_counts <- trt_est_from_sim %>% group_by(restype) %>% count(hasna = (TRUE==is.na(est_stderr)))
        num_calib_est_with_na_stderr <- na_stderr_counts %>%
            filter( (restype == "calibrated") & (hasna == TRUE)) %>% pull(n)

        sim_nums_w_na_calib_stderr <- trt_est_from_sim %>%
            filter(is.na(est_stderr) & restype == "calibrated") %>% pull(sim_num)

        trt_est_to_plot <- trt_est_from_sim %>% filter(!(sim_num %in% sim_nums_w_na_calib_stderr))

        sim_subtile <- glue(sim_desc, "\n",
                            "{length(sim_nums_w_na_calib_stderr)} studies removed due to ",
                            "NA in calibrated standard error")
    }

    # Calculates 95% CIs
    trt_est_from_sim_w_ci <- trt_est_to_plot %>%
        mutate(est_lower = est - (1.96 * est_stderr)) %>%
        mutate(est_upper = est + (1.96 * est_stderr))

    xlab_lower_limit <- quantile(trt_est_from_sim_w_ci %>% pull(est_lower), probs = 1-calibrated_max_est_se_quantile)[[1]]
    xlab_upper_limit <- quantile(trt_est_from_sim_w_ci %>% pull(est_upper), probs = calibrated_max_est_se_quantile)[[1]]

    sim_subtile <- glue(sim_subtile, "\n",
                        "Zooming in {calibrated_max_est_se_quantile} quantiles of x-axis")

    # Gets the true value
    true_eff_size <- trt_est_from_sim_w_ci %>% pull(true_eff) %>% unique()

    est_forest_plot <- trt_est_from_sim_w_ci %>% arrange(restype) %>%
        ggplot(aes(y = seq(1L, to = nrow(.)), x = est)) +
        geom_pointrange(aes(xmin = est_lower, xmax = est_upper, colour = restype), size = 0.2) +
        scale_colour_discrete(
            name = "Type",
            breaks = c("calibrated", "uncalibrated"),
            labels = c("Calibrated", "Uncalibrated")) +
        geom_vline(xintercept = true_eff_size) +
        ggtitle(plot_title, subtitle = sim_subtile) +
        ylab("Simulations") + xlab("Estimate")

    est_forest_plot <- est_forest_plot + coord_cartesian(xlim = c(xlab_lower_limit, xlab_upper_limit))

    return(est_forest_plot)

} # End plot_est_forest_plot()

#' Plots a general \emph{forest plot} with 95% confidence interval.
#'
#' @param data A \code{data.frame} that contains the estimates
#'    and standard error.
plot_forest_general <- function(data = NULL,
                                estimate_colname = "estimate",
                                stderr_colname = "std.error",
                                lower95_ci_colname = "lower95",
                                upper95_ci_colname = "upper95",
                                true_est_val = NA_real_,
                                aes_arrange_var = NA_character_,
                                calibrated_max_est_se_quantile = 1.0,
                                plot_title = "Forest plot of estimate of treatment effect",
                                plot_subtitle = "") {

    if (!require(dplyr))   { library(dplyr) }
    if (!require(ggplot2)) { library(ggplot2) }
    if (!require(glue))    { library(glue) }


    plot_data <- data

    num_sim_with_na_stderr <- NA_integer_

    # Calculates 95% CIs
    #trt_est_from_sim_w_ci <- plot_data %>%
    #    mutate(!!lower95_ci_colname := .data[[estimate_colname]] - (1.96 * .data[[stderr_colname]])) %>%
    #    mutate(!!upper95_ci_colname := .data[[estimate_colname]] + (1.96 * .data[[stderr_colname]]))
    # est_forest_plot <- trt_est_from_sim_w_ci %>% ggplot(aes(y = seq.int(1L, to = nrow(.))))

    if (!is.na(aes_arrange_var)) {
        est_forest_plot <- plot_data %>% arrange(!!sym(aes_arrange_var)) %>%
            ggplot(aes(y = seq.int(1L, to = nrow(.)))) +
            geom_pointrange(aes(x = !!sym(estimate_colname),
                                xmin = !!sym(lower95_ci_colname), xmax = !!sym(upper95_ci_colname),
                            colour = !!sym(aes_arrange_var)), size = .2, alpha = 0.7)
    } else {
        est_forest_plot <- plot_data %>%
            ggplot(aes(y = seq.int(1L, to = nrow(.)))) +
            geom_pointrange(aes(x = !!sym(estimate_colname),
                                xmin = !!sym(lower95_ci_colname), xmax = !!sym(upper95_ci_colname)),
                            size = .2, alpha = 0.7)
    }

    est_forest_plot <- est_forest_plot +
        ggtitle(plot_title, subtitle = plot_subtitle) +
        scale_colour_discrete(
            name = "Type",
            breaks = c("calibrated", "uncalibrated"),
            labels = c("Calibrated", "Uncalibrated")) +
        ylab("Simulations") + xlab("Estimate")

    if (!is.na(true_est_val)) {
        est_forest_plot <- est_forest_plot + geom_vline(xintercept = true_est_val)
    }

    return(est_forest_plot)


} # End plot_forest_general()


plot_bias_forest_plot <- function(bias_calc_sims = NULL,
                                 plot_title = "Forest plot of biases of treatment effect",
                                 sim_desc = "sim-desc",
                                 remove_sim_with_na_calib_stderr = FALSE) {

    trt_est_from_sim <- bias_calc_sims %>% tidyr::unnest(cols = c(res))

    num_calib_est_with_na_stderr <- NA_integer_

    trt_est_to_plot <- trt_est_from_sim
    sim_subtile <- sim_desc

    if (remove_sim_with_na_calib_stderr) {
        na_stderr_counts <- trt_est_from_sim %>% group_by(restype) %>% count(hasna = (TRUE==is.na(est_stderr)))
        num_calib_est_with_na_stderr <- na_stderr_counts %>%
            filter( (restype == "calibrated") & (hasna == TRUE)) %>% pull(n)

        sim_nums_w_na_calib_stderr <- trt_est_from_sim %>%
            filter(is.na(est_stderr) & restype == "calibrated") %>% pull(sim_num)

        trt_est_to_plot <- trt_est_from_sim %>% filter(!(sim_num %in% sim_nums_w_na_calib_stderr))

        sim_subtile <- glue(sim_desc, "\n",
                            "{length(sim_nums_w_na_calib_stderr)} studies removed due to ",
                            "NA in calibrated standard error")
    }

    bias_forest_plot <- trt_est_to_plot %>% arrange(restype) %>%
        ggplot(aes(y = seq(1L, to = nrow(.)), x = est)) +
        geom_point(aes(colour = restype)) +
        scale_colour_discrete(
            name = "Type",
            breaks = c("calibrated", "uncalibrated"),
            labels = c("Calibrated", "Uncalibrated")) +
        geom_vline(xintercept = 0.0) +
        ggtitle(plot_title, subtitle = sim_subtile) +
        ylab("") + xlab("Bias")

    return(bias_forest_plot)

} # End plot_bias_forest_plot()