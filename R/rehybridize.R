#' Virtual rehybridization with an array-to-sample assignment
#'
#' Create simulated dataset through "virtual rehybridization" for a given array-to-sample assignment.
#'
#' @param smp.eff the estimated sample effect dataset. The dataset must have rows as probes and columns as samples.
#' @param ary.eff the estimated array effect dataset. The dataset must have rows as probes and columns as samples.
#' It must have the same dimensions and the same probe names as the estimated sample effect dataset.
#' @param group.id a vector of sample-group labels for each sample of the estimated sample effect dataset.
#' It must be a 2-level non-numeric factor vector.
#' @param group.id.level a vector of sample-group label level. It must have two and only two elements and
#' the first element is the reference.
#' By default, \code{group.id.level = c("E", "V")}. That is in our study, we compare endometrial tumor samples to
#' ovarian tumor samples, with endometrial as our reference.
#' @param ary.to.smp.assign a vector of indices that assign arrays to samples (see details in \code{blocking.design},
#' \code{confounding.design} or \code{stratification.design}). It must have an equal length to the number of samples
#' in the estimated sample effect dataset.
#' The first half arrays in the vector have to be assigned to the sample group 1 and the second half to sample group 2.
#' @param icombat an indicator for combat adjustment. By default, \code{icombat = FALSE} for no ComBat adjustment.
#' @param isva an indicator for sva adjustment. By default, \code{isva = FALSE} for no sva adjustment.
#' @param iruv an indicator for RUV-4 adjustment. By default, \code{iruv = FALSE} for no RUV-4 adjustment.
#' @param smp.eff.ctrl the negative-control probe sample effect data if \code{iruv = TRUE}.
#' This dataset must have rows as probes and columns as samples.
#' It also must have the same number of samples and the same sample names as \code{smp.eff}.
#' @param ary.eff.ctrl the negative-control probe array effect data if \code{iruv = TRUE}.
#' It also must have the same dimensions and the same probe names as \code{smp.eff.ctrl}.
#' @return simulated data, after batch adjustment if specified
#' @import ruv sva
#' @importFrom stats model.matrix
#' @export
#' @keywords data.setup
#' @examples
#' \dontrun{
#' smp.eff <- estimate.smp.eff(uhdata = uhdata.pl)
#' ary.eff <- estimate.ary.eff(uhdata = uhdata.pl,
#'                              nuhdata = nuhdata.pl)
#'
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#'
#' smp.eff.nc <- smp.eff[!rownames(smp.eff) %in% ctrl.genes, ]
#' ary.eff.nc <- ary.eff[!rownames(ary.eff) %in% ctrl.genes, ]
#'
#' assign.ind <- confounding.design(seed = 1, num.smp = 192,
#' degree = "complete", rev.order = FALSE)
#'
#' group.id <- substr(colnames(smp.eff.nc), 7, 7)
#'
#' # no batch effect adjustment (default)
#' sim.data.raw <- rehybridize(smp.eff = smp.eff.nc,
#'                             ary.eff = ary.eff.nc,
#'                             group.id = group.id,
#'                             ary.to.smp.assign = assign.ind)
#'
#' # batch effect adjusting with sva
#' sim.data.sva <- rehybridize(smp.eff = smp.eff.nc,
#'                             ary.eff = ary.eff.nc,
#'                             group.id = group.id,
#'                             ary.to.smp.assign = assign.ind,
#'                             isva = TRUE)
#'
#' # batch effect adjusting with RUV-4
#' smp.eff.ctrl <- smp.eff[rownames(smp.eff) %in% ctrl.genes, ]
#' ary.eff.ctrl <- ary.eff[rownames(ary.eff) %in% ctrl.genes, ]
#'
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

