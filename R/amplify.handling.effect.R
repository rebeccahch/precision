#' Handling effect amplification
#'
#' Amplify handling effect in pre-specified slides by either a location shift or a scale change.
#'
#' @param handling.effect the estimated handling effect dataset to be modified. The dataset must have rows as probes and columns as samples.
#' @param amplify.array.id the array IDs specified to have its handling effect amplified.
#' If \code{type = "shift"} or \code{type = "scale1"}, a vector of array IDs must be supplied.
#' If \code{type = "scale2"}, a list of vectors of array IDs must be supplied;
#' each element in the list must be a vector of array IDs.
#' @param amplify.level a multiplier specified to amplify handling effect by.
#' A numeric multiplier must be supplied if \code{type = "shift"} or \code{type = "scale1"}.
#' A vector of multipliers must be supplied if type = "scale2" and it must have an equal length to the \code{amplify.array.id} list.
#' @param type a choice of amplification type, either "shift", "scale1" or "scale2" for either location shift
#' or scale change. By default, \code{type = "shift"}.
#' Location shift moves the entire specified arrays up or down by a constant.
#' Scale change 1 re-scales expressions that are in inter-quartiles towards the first and the third quartiles Within each array;
#' expressions that are outside of the inter-quartile range remain unchanged.
#' Scale change 2 re-scales expressions by the power of constants that are specified by the user for each batch.
#' @return an handling-effect-amplified set of handling effects
#' @keywords data.setup
#' @importFrom stats median
#' @export
#' @examples
#' \dontrun{
#' biological.effect <- estimate.biological.effect(uhdata = uhdata.pl)
#' handling.effect <- estimate.handling.effect(uhdata = uhdata.pl,
#'                             nuhdata = nuhdata.pl)
#'
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#'
#' biological.effect.nc <- biological.effect[!rownames(biological.effect) %in% ctrl.genes, ]
#' handling.effect.nc <- handling.effect[!rownames(handling.effect) %in% ctrl.genes, ]
#'
#' handling.effect.nc.tr <- handling.effect.nc[, c(1:64, 129:192)]
#'
#' # location shift
#' handling.effect.nc.tr.shift <- amplify.handling.effect(handling.effect = handling.effect.nc.tr,
#'                                        amplify.array.id = colnames(handling.effect.nc.tr)[1:64],
#'                                        amplify.level = 2, type = "shift")
#'
#' # scale change 1
#' handling.effect.nc.tr.scale1 <- amplify.handling.effect(handling.effect = handling.effect.nc.tr,
#'                                         amplify.array.id = colnames(handling.effect.nc.tr)[1:64],
#'                                         amplify.level = 2, type = "scale1")
#'
#' # scale change 2
#' amplify.array.id <- list(1:40, 41:64, (129:160) - 64, (161:192) - 64)
#' for(i in 1:length(amplify.array.id))
#'   amplify.array.id[[i]] <- colnames(handling.effect.nc.tr)[amplify.array.id[[i]]]
#' amplify.level <- c(1.2, 1.3, 1/3, 2/3)
#'
#' handling.effect.nc.tr.scale2 <- amplify.handling.effect(handling.effect = handling.effect.nc.tr,
#'                                         amplify.array.id = amplify.array.id,
#'                                         amplify.level = amplify.level,
#'                                         type = "scale2")
#'
#'
#' par(mfrow = c(2, 2), mar = c(4, 3, 2, 2))
#' rng <- range(handling.effect.nc.tr, handling.effect.nc.tr.shift,
#'              handling.effect.nc.tr.scale1, handling.effect.nc.tr.scale2)
#' boxplot(handling.effect.nc.tr, main = "original",
#'         ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
#' boxplot(handling.effect.nc.tr.shift, main = "shifted",
#'         ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
#' boxplot(handling.effect.nc.tr.scale1, main = "scaled 1",
#'         ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
#' boxplot(handling.effect.nc.tr.scale2, main = "scaled 2",
#'         ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
#' }
#'

"amplify.handling.effect" <- function(handling.effect,
                              amplify.array.id,
                              amplify.level,
                              type = "shift"){

  stopifnot(type %in% c("shift", "scale1", "scale2"))
  stopifnot(unlist(amplify.array.id) %in% colnames(handling.effect))
  stopifnot(is.numeric(unlist(amplify.level)))
  if(type == "scale2") stopifnot(unlist(amplify.array.id) %in% colnames(handling.effect))
  if(type == "scale2") stopifnot(length(amplify.level) == length(amplify.array.id))


  if(type == "shift"){
    a.e <- cbind(handling.effect[, colnames(handling.effect) %in% amplify.array.id] + amplify.level,
                    handling.effect[, !colnames(handling.effect) %in% amplify.array.id])
  } else if(type == "scale1") { # scale 1

    handling.effect.scaled <- handling.effect
    for(i in 1:length(amplify.array.id)) {
      slide <- handling.effect[, colnames(handling.effect) %in% amplify.array.id[i]]
      med <- median(slide)
      min <- min(slide)
      max <- max(slide)
      scaled <- ifelse(slide == med, slide,
                       ifelse(slide > med, slide + (max - slide)/amplify.level,
                              slide + (min - slide)/amplify.level))
      handling.effect.scaled[, colnames(handling.effect) %in% amplify.array.id[i]] <-
        scaled
    }

    a.e <- handling.effect.scaled

  } else { # scale 2
      handling.effect.scaled <- NULL
      for(i in 1:length(amplify.level)){
        x <- handling.effect[, colnames(handling.effect) %in% amplify.array.id[[i]]]
        handling.effect.scaled <- cbind(handling.effect.scaled,
                                   sign(x)*abs(x)^amplify.level[i])
      }
      handling.effect.scaled <- handling.effect.scaled[, colnames(handling.effect)]
      a.e <- handling.effect.scaled
    }

  return(a.e)
}

