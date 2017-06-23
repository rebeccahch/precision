#' Extracting errors from PRECISION (both non-FLEX and FLEX) output
#'
#' Extract all errors in the \code{error_store} list from PRECISION output produced from either \code{precision.simulate.flex}; reduce the PRECISION output for plotting.
#'
#' @param precision.obj.path the path and document name for the PRECISION object. By default, \code{precision.obj.path = NULL} indicates that the PRECISION object is already existed in the current global environment.
#' @param precision.obj.name a string for the name of the PRECISION R object.
#' @return a list of lists containing all the errors from the PRECISION object with the same names as the object
#' @keywords result.manage
#' @export
#' @examples
#' \dontrun{
#' ## Load + extract errors from a PRECISION object
#' precision.error.df <- extract.precision.error(precision.obj.path =
#'   "mypath/precision.result.Rdata",
#'   precision.obj.name = "precision.result")
#'   
#' ## Load + extract errors from a PRECISION FLEX object
#' precision.flex.error.df <- extract.precision.error(precision.obj.path =
#'   "mypath/precision.flex.result.Rdata",
#'   precision.obj.name = "precision.flex.result")
#'
#' ## Extract errors from a PRECISION FLEX object previously existed in global env
#' load("mypath/precision.a.Rdata") # load object named "precision.a" in global env
#' precision.flex.error.df <- extract.precision.error(precision.obj.name = "precision.flex.result")
#' }
#'
"extract.precision.error" <- function(precision.obj.path = NULL,
                                      precision.obj.name){

  if(!is.null(precision.obj.path)){
    load(precision.obj.path, envir = temp.env <- new.env())
    eval(parse(text = paste0("data <- temp.env$", precision.obj.name)))
  } else{
    eval(parse(text = paste0("data <- ", precision.obj.name)))
  }

  dat <- rep(list(NA), nrow(data$error_store))
  names(dat) <- rownames(data$error_store)
  dat <- apply(data$error_store, 1, function(x) x)
  # rename so that CC+ and CC- being the same "CC." after data.frame()
  for(i in 1:length(dat)) names(dat[[i]]) <- gsub("[+]", "1", 
                                                  gsub("[-]", "2", names(dat[[i]])))
  dat <- lapply(dat, data.frame)
  # rm(temp.env, data)
  return(dat)
}


