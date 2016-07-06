#' Array effect amplification
#'
#' Amplify array effect in pre-specified slides by either a location shift or a scale change.
#'
#' @param ary.eff the estimated array effect dataset to be modified. The dataset must have rows as probes and columns as samples.
#' @param amplify.ary.id the array IDs specified to have its array effect amplified.
#' If \code{type = "shift"} or \code{type = "scale1"}, a vector of array IDs must be supplied.
#' If \code{type = "scale2"}, a list of vectors of array IDs must be supplied.
#' @param amplify.level a multiplier specified to amplify array effect by.
#' A numeric multiplier must be supplied if \code{type = "shift"} or \code{type = "scale1"}.
#' A vector of multipliers must be supplied if type = "scale2" and it must have an equal length to the \code{amplify.ary.id} list.
#' @param type a choice of amplification type, either "shift", "scale1" or "scale2" for either location shift
#' or scale change. By default \code{type = "shift"}.
#' Location shift moves the entire specified arrays up or down by a constant.
#' Scale change 1 within each array, re-scales expressions that are in inter-quartiles towards the first and the third quartiles;
#' expressions that are outside of the inter-quartile range remain unchanged.
#' Scale change 2 re-scales the expressions by the power of constants that are specified by the user for each batch.
#' @return an array-effect-amplified set of array effects
#' @keywords data.setup
#' @importFrom stats median
#' @export
#' @examples
#' \dontrun{
#' smp.eff <- estimate.smp.eff(uhdata = uhdata.pl)
#' ary.eff <- estimate.ary.eff(uhdata = uhdata.pl,
#'                             nuhdata = nuhdata.pl)
#'
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#'
#' smp.eff.nc <- smp.eff[!rownames(smp.eff) %in% ctrl.genes, ]
#' ary.eff.nc <- ary.eff[!rownames(ary.eff) %in% ctrl.genes, ]
#'
#' ary.eff.nc.tr <- ary.eff.nc[, c(1:64, 129:192)]
#'
#' # location shift
#' ary.eff.nc.tr.shift <- amplify.ary.eff(ary.eff = ary.eff.nc.tr,
#'                                        amplify.ary.id = colnames(ary.eff.nc.tr)[1:64],
#'                                        amplify.level = 2, type = "shift")
#'
#' # scale change 1
#' ary.eff.nc.tr.scale1 <- amplify.ary.eff(ary.eff = ary.eff.nc.tr,
#'                                         amplify.ary.id = colnames(ary.eff.nc.tr)[1:64],
#'                                         amplify.level = 2, type = "scale1")
#'
#' # scale change 2
#' amplify.ary.id <- list(1:40, 41:64, (129:160) - 64, (161:192) - 64)
#' for(i in 1:length(amplify.ary.id))
#'   amplify.ary.id[[i]] <- colnames(ary.eff.nc.tr)[amplify.ary.id[[i]]]
#' amplify.level <- c(1.2, 1.3, 1/3, 2/3)
#'
#' ary.eff.nc.tr.scale2 <- amplify.ary.eff(ary.eff = ary.eff.nc.tr,
#'                                         amplify.ary.id = amplify.ary.id,
#'                                         amplify.level = amplify.level,
#'                                         type = "scale2")
#'
#'
#' par(mfrow = c(2, 2), mar = c(4, 3, 2, 2))
#' rng <- range(ary.eff.nc.tr, ary.eff.nc.tr.shift,
#'              ary.eff.nc.tr.scale1, ary.eff.nc.tr.scale2)
#' boxplot(ary.eff.nc.tr, main = "original",
#'         ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
#' boxplot(ary.eff.nc.tr.shift, main = "shifted",
#'         ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
#' boxplot(ary.eff.nc.tr.scale1, main = "scaled 1",
#'         ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
#' boxplot(ary.eff.nc.tr.scale2, main = "scaled 2",
#'         ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
#' }
#'

"amplify.ary.eff" <- function(ary.eff,
                              amplify.ary.id,
                              amplify.level,
                              type = "shift"){

  stopifnot(type %in% c("shift", "scale1", "scale2"))
  stopifnot(unlist(amplify.ary.id) %in% colnames(ary.eff))
  stopifnot(is.numeric(unlist(amplify.level)))
  if(type == "scale2") stopifnot(unlist(amplify.ary.id) %in% colnames(ary.eff))
  if(type == "scale2") stopifnot(length(amplify.level) == length(amplify.ary.id))


  if(type == "shift"){
    a.e <- cbind(ary.eff[, colnames(ary.eff) %in% amplify.ary.id] + amplify.level,
                    ary.eff[, !colnames(ary.eff) %in% amplify.ary.id])
  } else if(type == "scale1") { # scale 1

    ary.eff.scaled <- ary.eff
    for(i in 1:length(amplify.ary.id)) {
      slide <- ary.eff[, colnames(ary.eff) %in% amplify.ary.id[i]]
      med <- median(slide)
      min <- min(slide)
      max <- max(slide)
      scaled <- ifelse(slide == med, slide,
                       ifelse(slide > med, slide + (max - slide)/amplify.level,
                              slide + (min - slide)/amplify.level))
      ary.eff.scaled[, colnames(ary.eff) %in% amplify.ary.id[i]] <-
        scaled
    }

    a.e <- ary.eff.scaled

  } else { # scale 2
      ary.eff.scaled <- NULL
      for(i in 1:length(amplify.level)){
        x <- ary.eff[, colnames(ary.eff) %in% amplify.ary.id[[i]]]
        ary.eff.scaled <- cbind(ary.eff.scaled,
                                   sign(x)*abs(x)^amplify.level[i])
      }
      ary.eff.scaled <- ary.eff.scaled[, colnames(ary.eff)]
      a.e <- ary.eff.scaled
    }

  return(a.e)
}

