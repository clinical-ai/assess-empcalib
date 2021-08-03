# Calibrate controls using systematic error model.
#
# @param est_from_all_ctrls_tbl \code{data.frame} (or similar) data structure that contains
#     data on true value of estimate, estimated effect, and standard error of the estimated effect.
# @param sys_err_model Systematic error model from Empirical Calibration
# @param calib_est_w_true_eff_size The true effect size of data to calibrate. For negative controls
#     this would be value of zero.
# @param all_ctrls_tbl_coldef List of column definitions for \code{est_from_all_ctrls_tbl}
#     data structure.
# calibrate_control <- function(calib_est_w_true_eff_size = 0.0,
calibrate_control <- function(est_from_all_ctrls_tbl = NULL, sys_err_model = NULL,
                              all_ctrls_tbl_coldef = list(est_colname = "actual_eff_est",
                                                          se_est_colname = "actual_eff_est_stderr")) {

    if (!require(EmpiricalCalibration)) { library(EmpiricalCalibration) }
    if (!require(glue))                 { library(glue)  }
    if (!require(purrr))                { library(purrr) }


    if (is.null(est_from_all_ctrls_tbl) || is.null(sys_err_model)) {
        stop(glue("!! ERROR calibrate_control() NULL arguments received."))
    }

    ests_to_calibrate <- est_from_all_ctrls_tbl

    calibrate_est <- apply(ests_to_calibrate, MARGIN = 1, FUN = function(row) {

        estimate_to_calib <- row[[ all_ctrls_tbl_coldef[['est_colname']] ]]
        se_estimate_to_calib <- row[[ all_ctrls_tbl_coldef[['se_est_colname']] ]]

        calibrated_est_obj <- EmpiricalCalibration::calibrateConfidenceInterval(
            logRr = estimate_to_calib,
            seLogRr = se_estimate_to_calib,
            model = sys_err_model)

        return(c(actual_eff_est = calibrated_est_obj[['logRr']],
                 est_L95   = calibrated_est_obj[['logLb95Rr']],
                 est_U95   = calibrated_est_obj[['logUb95Rr']] ,
                 actual_eff_est_stderr     = calibrated_est_obj[['seLogRr']]))
    })

    calibrated_est_tbl <- calibrate_est %>% t() %>% as.data.frame()

    return(calibrated_est_tbl)

} # End calibrate_control()
