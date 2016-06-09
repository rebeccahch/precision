#' Prediction with nearest shrunken centroid classifier
#'
#' Predicts from a nearest shrunken centroid fit.
#'
#' @param pam.intcv.model a PAM classifier built with pam.intcv().
#' @param pred.obj expression dataset to have its sample group predicted, rows as probes, columns as samples; should have equal number of probes as the data trained.
#' @param pred.obj.group.id sample group corresponding to the data to be predicted; should have equal length as the number of samples as pred.obj.
#' @return predicted object, predicted error and predicted features
#' @export
#' @import pamr
#' @keywords classification
#' @examples
#' set.seed(101)
#' smp.eff <- estimate.smp.eff(r.data = r.data.pl)
#' ctrl.genes <- unique(rownames(r.data.pl))[grep("NC", unique(rownames(r.data.pl)))]
#' smp.eff.nc <- smp.eff[!rownames(smp.eff) %in% ctrl.genes, ]
#' group.id <- substr(colnames(smp.eff.nc), 7, 7)
#'
#' smp.eff.train.ind <- colnames(smp.eff.nc)[c(sample(which(group.id == "E"), size = 64),
#'                                           sample(which(group.id == "V"), size = 64))]
#' smp.eff.test.ind <- colnames(smp.eff.nc)[!colnames(smp.eff.nc) %in% smp.eff.train.ind]
#' smp.eff.nc.tr <- smp.eff.nc[, smp.eff.train.ind]
#' smp.eff.nc.te <- smp.eff.nc[, smp.eff.test.ind]
#'
#' pam.int <- pam.intcv(X = smp.eff.nc.tr,
#'                      y = substr(colnames(smp.eff.nc.tr), 7, 7),
#'                      kfold = 5, seed = 1)
#'
#' pam.pred <- pam.predict(pam.intcv.model = pam.int,
#'                         pred.obj = smp.eff.nc.te,
#'                         pred.obj.group.id = substr(colnames(smp.eff.nc.te), 7, 7))
#' pam.int$mc
#' pam.pred$mc

"pam.predict" <- function(pam.intcv.model, pred.obj, pred.obj.group.id){
  pred <- pamr::pamr.predict(pam.intcv.model$model, newx = pred.obj,
                       threshold = pam.intcv.model$model$threshold,
                       type = "class")

  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
  feature <- pamr::pamr.predict(pam.intcv.model$model, newx = pred.obj,
                                threshold = pam.intcv.model$model$threshold,
                                type = "posterior")
  return(list(pred=pred, mc=mc, feature=feature))
}