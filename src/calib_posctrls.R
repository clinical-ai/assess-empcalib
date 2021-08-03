calib_posctrl <- function(ctrl_est_tbl = NULL, err_model = NULL,
                          use_adj_true_colname = FALSE,
                          ctrl_est_tbl_cols_def = list(sim_id_colname = "sim_id",
                                                       outcome_id_colname = "outcome_id",
                                                       negctrl_id_colname = "negctrl_id",
                                                       true_eff_colname = "true_eff",
                                                       adj_true_colname = "adj_true_eff",
                                                       est_colname = "actual_eff_est",
                                                       se_est_colname = "actual_eff_est_stderr")) {

    ctrl_data_to_calibrate <- ctrl_est_tbl

    calibrated_posctrl_data <- calibrate_control(
        est_from_all_ctrls_tbl = ctrl_data_to_calibrate,
        sys_err_model = err_model,
        all_ctrls_tbl_coldef = list(est_colname = ctrl_est_tbl_cols_def[['est_colname']],
                                    se_est_colname = ctrl_est_tbl_cols_def[['se_est_colname']] ))

    ctrl_est_data_for_combine <- ctrl_data_to_calibrate %>% select(c(
        ctrl_est_tbl_cols_def[['sim_id_colname']],
        ctrl_est_tbl_cols_def[['outcome_id_colname']],
        ctrl_est_tbl_cols_def[['negctrl_id_colname']],
        ctrl_est_tbl_cols_def[['true_eff_colname']] ))

    if (use_adj_true_colname) {
        ctrl_est_data_for_combine <- ctrl_est_data_for_combine %>% mutate(
        !!sym(ctrl_est_tbl_cols_def[['true_eff_colname']]) := ctrl_est_tbl %>% pull(ctrl_est_tbl_cols_def[['adj_true_colname']]) )
    }

    calibrated_result <- cbind(ctrl_est_data_for_combine , calibrated_posctrl_data) %>%
        mutate(restype = rep("calibrated", times = nrow(calibrated_posctrl_data)))

    if (use_adj_true_colname) {
        ctrl_data_to_calibrate <- ctrl_data_to_calibrate %>% mutate(
            !!sym(ctrl_est_tbl_cols_def[['true_eff_colname']]) := ctrl_est_tbl %>% pull(ctrl_est_tbl_cols_def[['adj_true_colname']]) )
    }

    uncalibrated_ctrl_data <- ctrl_data_to_calibrate %>%
        mutate(restype = rep("uncalibrated", times = nrow(ctrl_data_to_calibrate)))

    tbl_to_rt <- merge_calib_uncalib_ctrl_data(uncalibrated_ctrl_data, calibrated_result)

} # End defining calib_posctrl()
