#' Estimated Sample Effects
#'
#' Estimate biological effects of a sample from the uniformly-handled dataset.
#'
#' @param uhdata the uniformly-handled dataset.
#' The dataset must have rows as probes and columns as samples.
#' @return an estimation of the sample effects
#' @keywords data.setup
#' @export
#' @examples
#' biological.effect <- estimate.biological.effect(uhdata = uhdata.pl)
#'

"estimate.biological.effect" <- function(uhdata) {
    biological.effect <- uhdata
    return(biological.effect)
  }
