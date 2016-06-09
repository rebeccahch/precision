#' Rehybridization with an array-to-sample assignment
#'
#' Creates simulated dataset through rehybridization with a specified array-to-sample assignment.
#'
#' @param smp.eff sample effect data, rows as probes, columns as samples.
#' @param ary.eff array effect data, rows as probes, columns as samples; must have same dimensions and same probe name as sample effect data.
#' @param group.id sample group label; must be a 2-level non-numeric factor vector.
#' @param group.id.level sample group label level, the first one being the reference level; default = c("E", "V") in our studies when comparing endometrial to ovarian samples.
#' @param ary.to.smp.assign array-to-sample assignment, equal length as number of samples of sample effect data; first half of the vector assigning to endometrial, second half to ovarian.
#' @param icombat indicator for combat adjustment; default is not to adjust, icombat = FALSE.
#' @param isva indicator for sva adjustment; default is not to adjust, isva = FALSE.
#' @param iruv indicator for RUV-4 adjustment; default is not to adjust, iruv = FALSE.
#' @param smp.eff.ctrl negative-control gene sample effect data if iruv = TRUE.
#' @param ary.eff.ctrl negative-control gene array effect data if iruv = TRUE
#' @return simulated data, after batch adjustment if specified
#' @import ruv
#' @import sva
#' @importFrom stats model.matrix
#' @export
#' @keywords data.setup
#' @examples
#' \dontrun{
#' smp.eff <- estimate.smp.eff(r.data = r.data.pl)
#' ary.eff <- estimate.ary.eff(r.data = r.data.pl,
#'                              non.r.data = non.r.data.pl)
#' ctrl.genes <- unique(rownames(r.data.pl))[grep("NC", unique(rownames(r.data.pl)))]
#' smp.eff.nc <- smp.eff[!rownames(smp.eff) %in% ctrl.genes, ]
#' ary.eff.nc <- ary.eff[!rownames(ary.eff) %in% ctrl.genes, ]
#' assign.ind <- confounding.design(seed = 1, num.smp = 192,
#' degree = "complete", rev.order = FALSE)
#' group.id <- substr(colnames(smp.eff.nc), 7, 7)
#' sim.data.raw <- rehybridize(smp.eff = smp.eff.nc,
#'                             ary.eff = ary.eff.nc,
#'                             group.id = group.id,
#'                             ary.to.smp.assign = assign.ind)
#' sim.data.sva <- rehybridize(smp.eff = smp.eff.nc,
#'                             ary.eff = ary.eff.nc,
#'                             group.id = group.id,
#'                             ary.to.smp.assign = assign.ind,
#'                             isva = TRUE)
#' smp.eff.ctrl <- smp.eff[rownames(smp.eff) %in% ctrl.genes, ]
#' ary.eff.ctrl <- ary.eff[rownames(ary.eff) %in% ctrl.genes, ]
#' sim.data.ruv <- rehybridize(smp.eff = smp.eff.nc,
#'                             ary.eff = ary.eff.nc,
#'                             group.id = group.id,
#'                             ary.to.smp.assign = assign.ind,
#'                             iruv = TRUE,
#'                             smp.eff.ctrl = smp.eff.ctrl,
#'                             ary.eff.ctrl = ary.eff.ctrl)
#' }

"rehybridize" <- function (smp.eff,
                                ary.eff,
                                group.id,
                                group.id.level = c("E", "V"),
                                ary.to.smp.assign,
                                icombat = FALSE,
                                isva = FALSE,
                                iruv = FALSE,
                                smp.eff.ctrl = NULL,
                                ary.eff.ctrl = NULL) {
  stopifnot(dim(smp.eff) == dim(ary.eff))
  stopifnot(rownames(smp.eff) == rownames(ary.eff))
  stopifnot(group.id %in% group.id.level)
  stopifnot(sum(icombat, isva, iruv) < 2)
  if(iruv) stopifnot(!is.null(smp.eff.ctrl) & !is.null(ary.eff.ctrl))

  halfcut <- length(ary.to.smp.assign)/2
  out <- cbind(smp.eff[, group.id == group.id.level[1]] +
                 ary.eff[, ary.to.smp.assign[1:halfcut]],
               smp.eff[, group.id == group.id.level[2]] +
                 ary.eff[, ary.to.smp.assign[(halfcut + 1):length(ary.to.smp.assign)]])

  if(icombat){ # ComBat
    cat("ComBat adjusting\n")
    mod.tr = model.matrix(~rep(1, ncol(out)))
    batch = cut(ary.to.smp.assign, 8*(0:(ncol(out)/8))) # adjust by array slide
    table(batch, substr(colnames(out), 7, 7))

    combat_dat = sva::ComBat(dat = out,
                        batch = batch,
                        mod = mod.tr,
                        par.prior = TRUE)
    out <- combat_dat
  }

  if(isva){ # sva
    mod0 = model.matrix(~1, data = data.frame(colnames(out)))
    mod = model.matrix(~rep(c(0, 1), each = halfcut)) # current "out" is not in original order
    n.sv = sva::num.sv(dat = out, mod = mod, method = "leek")
    sva_dat = sva::sva(out, mod, mod0, n.sv=n.sv)
    out <- out[, colnames(smp.eff)]
    return(list(trainData = out, trainMod = mod, trainSV = sva_dat))
    stop()
  }

  if(iruv){ # ruv
    halfcut <- length(ary.to.smp.assign)/2
    ctrl <- cbind(smp.eff.ctrl[, group.id == group.id.level[1]] +
                    ary.eff.ctrl[, ary.to.smp.assign[1:halfcut]],
                  smp.eff.ctrl[, group.id == group.id.level[2]] +
                    ary.eff.ctrl[, ary.to.smp.assign[(halfcut + 1):length(ary.to.smp.assign)]])
    combine_dat <- rbind(out, ctrl)
    ctrl.ind <- rownames(combine_dat) %in% unique(rownames(ctrl))

    cat("RUV4 normalize \n")
    Y <- t(combine_dat)
    X <- rep(c(0, 1), each = halfcut) # current "combine_dat" is not in original order
    temp <- ruv::getK(Y = Y, X = matrix(X), ctl = ctrl.ind)
    xx <- ruv::RUV4(Y = Y, X = matrix(X), ctl = ctrl.ind, k = temp$k)
    ruv_dat <- Y - xx$W %*% solve(t(xx$W) %*% xx$W, t(xx$W) %*% Y) # no shrinkage
    rm(temp, xx, X, Y)
    ruv_dat <- t(ruv_dat)[!ctrl.ind, ]

    out <- ruv_dat
  }

  out <- out[, colnames(smp.eff)]
  return(out)
}

