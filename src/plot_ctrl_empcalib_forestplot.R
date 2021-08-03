#' Plots funnel plot of calibrated (and un-calibrated) control data.
#' Similar to
#' https://ohdsi.github.io/EmpiricalCalibration/articles/EmpiricalCiCalibrationVignette.html#confidence-interval-calibration-1
#' Using function interface in v2.1.0 of EmpiricalCalibration package.
#'
#' @param calib_data \code{\link{data.frame}} like structure containing the calibrated and
#'     un-calibrated estimates of a control.
#' @param calib_data_cols Named list with column names of relevant data in \code{calib_data}
#' @param syserr_model Systematic error model returned from \code{\link{EmpiricalCalibration::fitSystematicErrorModel}}.
#'     If \code{NULL}, \link{\code{EmpiricalCalibration::plotCiCalibrationEffect}} will attempt to
#'     fit one.
#' @param syserr_model_uselegacy Defaults to \code{FALSE}. Parameter when fitting a systematic
#'     error model in \link{\code{EmpiricalCalibration::plotCiCalibrationEffect}}.
#'
#' @return \link{\code{ggplot}} object.
#'
plot_control_empclib_funnelplot <- function(calib_data = NULL,
                                            calib_data_cols = list(
                                                est_colname = "actual_eff_est",
                                                est_stderr_colname = "actual_eff_est_stderr",
                                                true_eff_colname = "true_eff",
                                                adjusted_true_eff_colname = "adj_true_eff",
                                                calib_indicator_colname = "restype"
                                            ),
                                            syserr_model = NULL,
                                            syserr_model_uselegacy = FALSE,
                                            title = NULL,
                                            adj_true_eff_posctrl = FALSE) {

    if (!require(EmpiricalCalibration)) { library(EmpiricalCalibration) }

    true_eff_vec <- calib_data[[ calib_data_cols[['true_eff_colname']] ]]
    if (adj_true_eff_posctrl) {
        true_eff_vec <- calib_data[[ calib_data_cols[['adjusted_true_eff_colname']] ]]
    }
    calib_eff_funnel_plot <- plotCiCalibrationEffect(
        logRr     = calib_data[[ calib_data_cols[['est_colname']] ]],
        seLogRr   = calib_data[[ calib_data_cols[['est_stderr_colname']] ]],
        trueLogRr = true_eff_vec,
        legacy    = syserr_model_uselegacy,
        model     = syserr_model,
        xLabel    = "Log Odds Ratio",
        title     = title
    )

    return(calib_eff_funnel_plot)

} # End plot_control_empclib_funnelplot()

plot_calib_n_uncalib_control_empclib_funnelplots <- function(calib_data = NULL,
                                                             calib_data_cols = list(
                                                                 est_colname = "actual_eff_est",
                                                                 est_stderr_colname = "actual_eff_est_stderr",
                                                                 true_eff_colname = "true_eff",
                                                                 adjusted_true_eff_colname = "adj_true_eff",
                                                                 calib_indicator_colname = "restype",
                                                                 calib_val = "calibrated",
                                                                 uncalib_val = "uncalibrated"
                                                             ),
                                                             syserr_model = NULL,
                                                             syserr_model_uselegacy = FALSE,
                                                             title = NULL,
                                                             posctrl_data = FALSE,
                                                             nposctrls = NA_integer_,
                                                             adj_true_eff_posctrl = FALSE) {

    if (!require(dplyr)) { library(dplyr) }

    # Obtain funnel plot for un-calibrated data
    calib_data_for_plot <- calib_data %>%
        filter(!!sym(calib_data_cols[['calib_indicator_colname']]) == calib_data_cols[['uncalib_val']])

    uncalib_ctrl_funnel_plot <- plot_control_empclib_funnelplot(calib_data_for_plot,
                                                                syserr_model = syserr_model,
                                                                syserr_model_uselegacy = syserr_model_uselegacy,
                                                                title = title,
                                                                adj_true_eff_posctrl = adj_true_eff_posctrl)

    # Obtain funnel plot for calibrated data
    calib_data_for_plot <- calib_data %>%
        filter(!!sym(calib_data_cols[['calib_indicator_colname']]) == calib_data_cols[['calib_val']])
    calib_ctrl_funnel_plot <- plot_control_empclib_funnelplot(calib_data_for_plot,
                                                              syserr_model = syserr_model,
                                                              syserr_model_uselegacy = syserr_model_uselegacy,
                                                              title = title,
                                                              adj_true_eff_posctrl = adj_true_eff_posctrl)

    # Combines the plots
    if (!require(ggpubr)) { library(ggpubr) }

    ncol_combined_plot <- 1L # Default one, but for positive controls, let user specify
    if (posctrl_data) {
        ncol <- nposctrls
    }
    combined_funnel_plot <- ggarrange(uncalib_ctrl_funnel_plot, calib_ctrl_funnel_plot,
              labels = c("Uncalibrated", "Calibrated"),
              ncol = ncol_combined_plot, nrow = 2)

    return(combined_funnel_plot)

} # End plot_calib_n_uncalib_control_empclib_funnelplots()
