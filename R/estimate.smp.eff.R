#' Estimated Sample Effects
#'
#' Estimate sample effects from the expressions of the uniformly-handled data.
#'
#' @param uhdata the uniformly-handled expression dataset.
#' The dataset must have rows as probes and columns as samples.
#' @return an estimation of the sample effects
#' @keywords data.setup
#' @export
#' @examples
#' smp.eff <- estimate.smp.eff(uhdata = uhdata.pl)
#'

"estimate.smp.eff" <- function(uhdata) {
    smp.eff <- uhdata
    return(smp.eff)
  }
