load_config <- function(csvconfigfile = NULL,
                        configcols = list(setting_colname = "flag", value_colname = "value")) {

    if (!require(dplyr)) { library(dplyr) }
    if (!require(readr)) { library(readr) }
    if (!require(stringr)) { library(stringr) }

    if (is.null(csvconfigfile) | (0 == str_length(csvconfigfile))) { return (NULL) }

    user_settings <- readr::read_csv(csvconfigfile)
    setting_col <- configcols[['setting_colname']]
    val_col <- configcols[['value_colname']]

    return(function(setting_to_get = NULL, default_val = NULL) {
        param <- user_settings %>%
            filter(!!sym(setting_col) == setting_to_get) %>% pull(!!sym(val_col))

        param_to_return <- param
        # Empty param
        if (0 == length(param)) { param_to_return <- default_val }

        return(param_to_return)
    }) # End inner function
} # End load_config() function
