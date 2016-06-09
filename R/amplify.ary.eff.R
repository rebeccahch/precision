#' Array effect amplification
#'
#' Amplifies array effect in specified slides in a training set by a multiplier.
#'
#' @param ary.eff.tr the array effect training set to be modified, rows as probes, columns as samples.
#' @param amplify.slide.id the slide IDs specified to have its array effect amplified; a vector of slide IDs if type = "shift" or "scale1", a list of vectors of slide IDs if type = "scale2".
#' @param amplify.level a multiplier specified to amplify array effect by; a multiplier if type = "shift" or "scale1"; a vector of multipliers if type = "scale2" which has to be equaal length to the amplify.slide.id list.
#' @param type a choice of amplification type, either "shift", "scale1" or "scale2" for either location shift
#' amplification or scale amplification; default is "shift". Location shift amplification shifts the entire
#' specified arrays up by a constant. Scaling amplification scales the expressions towards maximum and minimum per array.
#' @return a array-effect-amplified array effect training data
#' @keywords data.setup
#' @importFrom stats median
#' @export
#' @examples
#' \dontrun{
#' smp.eff <- estimate.smp.eff(r.data = r.data.pl)
#' ary.eff <- estimate.ary.eff(r.data = r.data.pl,
#'                              non.r.data = non.r.data.pl)
#' ctrl.genes <- unique(rownames(r.data.pl))[grep("NC", unique(rownames(r.data.pl)))]
#' smp.eff.nc <- smp.eff[!rownames(smp.eff) %in% ctrl.genes, ]
#' ary.eff.nc <- ary.eff[!rownames(ary.eff) %in% ctrl.genes, ]
#' ary.eff.nc.tr <- ary.eff.nc[, c(1:64, 129:192)]
#' ary.eff.nc.tr.shift <- amplify.ary.eff(ary.eff.tr = ary.eff.nc.tr,
#'                                        amplify.slide.id = colnames(ary.eff.nc.tr)[1:64],
#'                                        amplify.level = 2, type = "shift")
#' ary.eff.nc.tr.scale1 <- amplify.ary.eff(ary.eff.tr = ary.eff.nc.tr,
#'                                         amplify.slide.id = colnames(ary.eff.nc.tr)[1:64],
#'                                         amplify.level = 2, type = "scale1")
#' amplify.slide.id <- list(1:40, 41:64, (129:160) - 64, (161:192) - 64)
#' for(i in 1:length(amplify.slide.id))
#'   amplify.slide.id[[i]] <- colnames(ary.eff.nc.tr)[amplify.slide.id[[i]]]
#' amplify.level <- c(1.2, 1.3, 1/3, 2/3)
#' ary.eff.nc.tr.scale2 <- amplify.ary.eff(ary.eff.tr = ary.eff.nc.tr,
#'                                         amplify.slide.id = amplify.slide.id,
#'                                         amplify.level = amplify.level,
#'                                         type = "scale2")
#' par(mfrow = c(2, 2), mar = c(4, 3, 2, 2))
#' rng <- range(ary.eff.nc.tr, ary.eff.nc.tr.shift, ary.eff.nc.tr.scale1, ary.eff.nc.tr.scale2)
#' boxplot(ary.eff.nc.tr, main = "original",
#'         ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
#' boxplot(ary.eff.nc.tr.shift, main = "shifted",
#'         ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
#' boxplot(ary.eff.nc.tr.scale1, main = "scaled 1",
#'         ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
#' boxplot(ary.eff.nc.tr.scale2, main = "scaled 2",
#'         ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
#'}


"amplify.ary.eff" <- function(ary.eff.tr,
                              amplify.slide.id,
                              amplify.level,
                              type = "shift"){

  stopifnot(type %in% c("shift", "scale1", "scale2"))
  stopifnot(unlist(amplify.slide.id) %in% colnames(ary.eff.tr))
  stopifnot(is.numeric(unlist(amplify.level)))
  if(type == "scale2") stopifnot(unlist(amplify.slide.id) %in% colnames(ary.eff.tr))
  if(type == "scale2") stopifnot(length(amplify.level) == length(amplify.slide.id))


  if(type == "shift"){
    a.e.tr <- cbind(ary.eff.tr[, colnames(ary.eff.tr) %in% amplify.slide.id] + amplify.level, ary.eff.tr[, !colnames(ary.eff.tr) %in% amplify.slide.id])
  } else if(type == "scale1") { # scale 1

    ary.eff.tr.scaled <- ary.eff.tr
    for(i in 1:length(amplify.slide.id)) {
      slide <- ary.eff.tr[, colnames(ary.eff.tr) %in% amplify.slide.id[i]]
      med <- median(slide)
      min <- min(slide)
      max <- max(slide)
      scaled <- ifelse(slide == med, slide,
                       ifelse(slide > med, slide + (max - slide)/amplify.level,
                              slide + (min - slide)/amplify.level))
      ary.eff.tr.scaled[, colnames(ary.eff.tr) %in% amplify.slide.id[i]] <-
        scaled
    }

    a.e.tr <- ary.eff.tr.scaled

  } else { # scale 2
      ary.eff.tr.scaled <- NULL
      for(i in 1:length(amplify.level)){
        x <- ary.eff.tr[, colnames(ary.eff.tr) %in% amplify.slide.id[[i]]]
        ary.eff.tr.scaled <- cbind(ary.eff.tr.scaled,
                                   sign(x)*abs(x)^amplify.level[i])
      }
      ary.eff.tr.scaled <- ary.eff.tr.scaled[, colnames(ary.eff.tr)]
      a.e.tr <- ary.eff.tr.scaled
    }

  return(a.e.tr)
}

