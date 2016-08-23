#' Estimated Sample Effects
#'
#' Estimate biological effects of a sample from the uniformly-handled dataset.
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
