#' Plots biases of calibrated and uncalibrated results.
#'
#' @param bias_calc_sims Table of biases obtained from calling
#'     \code{extract_tbl_from_sim_result()} function.
#'
plot_bias_calc_boxplot <- function(bias_calc_sims = NULL,
                                   plot_title = "Boxplot of biases",
                                   sim_desc = "sim-desc",
                                   axis_break_gap = 0.1) {

    if (!require(tidyr)) { library(tidyr) }
    if (!require(ggplot2)) { library(ggplot2) }

    bias_calc_data <- bias_calc_sims %>% tidyr::unnest(res)

    y_axis_breaks <- seq(floor(min(bias_calc_data$est_bias)),
                         ceiling(max(bias_calc_data$est_bias)),
                         by = axis_break_gap)
    bias_calc_plot <- bias_calc_data %>% ggplot(aes(x = restype, y = est_bias, fill = restype)) +
        geom_boxplot(width = 0.2) +
        scale_fill_discrete(
            name = "Type",
            breaks = c("calibrated", "uncalibrated"),
            labels = c("Calibrated", "Uncalibrated")) +
        scale_y_continuous(breaks = y_axis_breaks) +
        geom_hline(yintercept = 0.0, linetype = "dotted", color = "blue") +
        ggtitle(plot_title, subtitle = sim_desc) +
        xlab("Estimate type") +
        ylab("Bias (log odds ratio)")

    return(bias_calc_plot)

} # End plot_bias_calc()

# Zooms into the boxplot of biases, removing outlier estimates from view
zoom_in_bias_calc_boxplot <- function(bias_boxplot = NULL, bias_calc_sims = NULL,
                                      quantile_prob_limits = c(0.05, 0.95),
                                      bias_colname = "est_bias",
                                      plot_title = "Boxplot of biases zoomed 0.05-0.95 quantiles",
                                      sim_desc = "sim_desc") {

    calibrated_est_range_for_plotting <- bias_calc_sims %>% tidyr::unnest(res) %>%
        filter(restype == "calibrated") %>% pull(.data[[bias_colname]]) %>% quantile(prob = quantile_prob_limits)
    uncalibrated_est_range_for_plotting <- bias_calc_sims %>% tidyr::unnest(res) %>%
        filter(restype == "uncalibrated") %>% pull(.data[[bias_colname]]) %>% quantile(prob = quantile_prob_limits)

    est_range_for_plotting <- rbind(calibrated_est_range_for_plotting, uncalibrated_est_range_for_plotting)
    est_lower_bound_for_plotting <- est_range_for_plotting[, 1] %>% min()
    est_upper_bound_for_plotting <- est_range_for_plotting[, 2] %>% max()

    # Zoom into the box plot
    zoomed_boxplot <- bias_boxplot +
        coord_cartesian(ylim = c(est_lower_bound_for_plotting, est_upper_bound_for_plotting))
    zoomed_boxplot <- zoomed_boxplot + ggtitle(plot_title, subtitle = sim_desc)

} # End zoom_in_bias_calc_boxplot()


#' Plots \emph{histogram} biases of calibrated and uncalibrated results.
#'
#' @param bias_calc_sims Table of biases obtained from calling
#'     \code{extract_tbl_from_sim_result()} function.
#'
plot_bias_calc_histogram <- function(bias_calc_sims = NULL, sim_desc = "sim-desc") {

    if (!require(tidyr)) { library(tidyr) }
    if (!require(ggplot2)) { library(ggplot2) }

    bias_calc_data <- bias_calc_sims %>% tidyr::unnest(res)

    bias_histogram <- bias_calc_data %>% ggplot(aes(x = est_bias, fill = restype), binwidth = "0.70") +
        geom_histogram(alpha = 0.6, bins = 60) +
        scale_fill_discrete(name = "Type",
                            breaks = c("calibrated", "uncalibrated"),
                            labels = c("Calibrated", "Uncalibrated")) +
        ggtitle("Histogram of bias in estimate of treatment effect", subtitle = sim_desc)
    bias_histogram <- bias_histogram + geom_vline(xintercept = 0.0)

    # Vertical lines for mean of estimates
    mean_calibrated_bias <- bias_calc_data %>% filter(restype == "calibrated") %>%
        pull(est_bias) %>% mean()
    mean_uncalibrated_bias <- bias_calc_data %>% filter(restype == "uncalibrated") %>%
        pull(est_bias) %>% mean()

    est_means_df <- data.frame(restype = c("calibrated", "uncalibrated"),
                            bias = c(mean_calibrated_bias, mean_uncalibrated_bias))
    bias_histogram <- bias_histogram + geom_vline(data = est_means_df,
        aes(xintercept = bias, linetype = restype))
    bias_histogram <- bias_histogram + scale_linetype_manual(
        name = "Mean of biases",
        breaks = c("calibrated", "uncalibrated"),
        values = c("dashed", "dotted"),
        labels = c("mean of calibrated", "mean of uncalibrated"))

    #bias_histogram <- bias_histogram + geom_vline(xintercept = mean_uncalibrated_bias, linetype = "dashed")
    #bias_histogram <- bias_histogram + geom_vline(xintercept = mean_calibrated_bias, linetype = "dotted")


    return(bias_histogram)

} # End plot_bias_calc_histogram()

#' Plots \emph{density} of biases of calibrated and uncalibrated results.
#' Suitable for large simulations
#'
#' @param bias_calc_sims Table of biases obtained from calling
#'     \code{extract_tbl_from_sim_result()} function.
#'
plot_bias_calc_density <- function(bias_calc_sims = NULL, sim_desc = "sim-desc") {

    if (!require(tidyr)) { library(tidyr) }
    if (!require(ggplot2)) { library(ggplot2) }

    bias_calc_data <- bias_calc_sims %>% tidyr::unnest(res)

    bias_density_plot <- bias_calc_data %>% ggplot(aes(x = est_bias, fill = restype), binwidth = "0.70") +
        geom_density(alpha = 0.6) +
        geom_rug(aes(colour = restype), sides = "b") +
        scale_fill_discrete(
            name = "Type", breaks = c("calibrated", "uncalibrated"),
            labels = c("Calibrated", "Uncalibrated")) +
        scale_color_discrete(
            name = "Type", breaks = c("calibrated", "uncalibrated"),
            labels = c("Calibrated", "Uncalibrated")) +
        ggtitle("Density of bias in estimate of treatment effect", subtitle = sim_desc)
    bias_density_plot <- bias_density_plot + geom_vline(xintercept = 0.0)

    # Vertical lines for mean of estimates
    mean_calibrated_bias <- bias_calc_data %>% filter(restype == "calibrated") %>%
        pull(est_bias) %>% mean()
    mean_uncalibrated_bias <- bias_calc_data %>% filter(restype == "uncalibrated") %>%
        pull(est_bias) %>% mean()

    est_means_df <- data.frame(restype = c("calibrated", "uncalibrated"),
                               bias = c(mean_calibrated_bias, mean_uncalibrated_bias))
    bias_density_plot <- bias_density_plot + geom_vline(data = est_means_df,
                                                  aes(xintercept = bias, linetype = restype))
    bias_density_plot <- bias_density_plot + scale_linetype_manual(
        name = "Mean of biases",
        breaks = c("calibrated", "uncalibrated"),
        values = c("dashed", "dotted"),
        labels = c("Mean of calibrated", "Mean of uncalibrated"))

    #bias_density_plot <- bias_density_plot + geom_vline(xintercept = mean_uncalibrated_bias, linetype = "dashed")
    #bias_density_plot <- bias_density_plot + geom_vline(xintercept = mean_calibrated_bias, linetype = "dotted")


    return(bias_density_plot)

} # End plot_bias_calc_density()
