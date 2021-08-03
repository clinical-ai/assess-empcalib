## ----------------------------------------------------------------------------
##
## Functional programming
##
## ----------------------------------------------------------------------------

if_cond_do_general <- function(data = NULL, condition = FALSE, fun = NULL, fun_args = list()) {
    if (TRUE == condition) {
        do.call(what = fun, args = unlist(list(data, fun_args) , recursive = FALSE))
    } else {
        return(data)
    }
}

if_cond_do <- function(data = NULL, condition = FALSE, fun = NULL) {
    if (TRUE == condition) {
        do.call(what = fun, args = list(data))
    } else {
        return(data)
    }
}

add_fn <- function(val) { return(function(x) { x + val } ) }
