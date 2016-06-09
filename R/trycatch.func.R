"trycatch.func" <- function(expr, msg = "") {
  out <- tryCatch({
    expr

  }, warning = function(cond) {
    message("warning: ")
    message(cond)
    # Choose a return value in case of warning
    return(NULL)

  }, error = function(cond) {
    message("error: ")
    message(cond)
    # Choose a return value in case of error
    return(NA)

  }, finally={
    message(msg)

  })
  return(out)
}
