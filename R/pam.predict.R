#' Prediction with nearest shrunken centroid classifier
#'
#' Predict from a nearest shrunken centroid fit.
#'
#' @references T. Hastie, R. Tibshirani, Balasubramanian Narasimhan and Gil Chu (2014).
#' pamr: Pam: prediction analysis for microarrays. R package version 1.55.
#' https://CRAN.R-project.org/package=pamr
#' @param pam.intcv.model a PAM classifier built with \code{pam.intcv()}.
#' @param pred.obj expression dataset to have its sample group predicted.
#' The dataset must have rows as probes and columns as samples.
#' It must have an equal number of probes as the dataset being trained.
#' @param pred.obj.group.id a vector of sample-group labels for
#' each sample of the dataset to be predicted.
#' It must have an equal length to the number of samples as \code{pred.obj}.
#' @return a list of 3 elements:
#' \item{pred}{predicted sample group for each sample}
#' \item{mc}{a predicted misclassification error rate (external validation)}
#' \item{prob}{predicted probability for each sample}
#' @export
#' @import pamr
#' @keywords classification
#' @examples
#' set.seed(101)
#' smp.eff <- estimate.smp.eff(uhdata = uhdata.pl)
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#' smp.eff.nc <- smp.eff[!rownames(smp.eff) %in% ctrl.genes, ]
#' group.id <- substr(colnames(smp.eff.nc), 7, 7)
#'
#' smp.eff.train.ind <- colnames(smp.eff.nc)[c(sample(which(group.id == "E"), size = 64),
#'                                           sample(which(group.id == "V"), size = 64))]
#' smp.eff.test.ind <- colnames(smp.eff.nc)[!colnames(smp.eff.nc) %in% smp.eff.train.ind]
#'
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
#'

"pam.predict" <- function(pam.intcv.model, pred.obj, pred.obj.group.id){
  pred <- pamr::pamr.predict(pam.intcv.model$model, newx = pred.obj,
                       threshold = pam.intcv.model$model$threshold,
                       type = "class")

  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
  prob <- pamr::pamr.predict(pam.intcv.model$model, newx = pred.obj,
                                threshold = pam.intcv.model$model$threshold,
                                type = "posterior")
  return(list(pred=pred, mc=mc, prob=prob))
}
