#' Prediction with nearest shrunken centroid classifier
#'
#' Predict from a nearest shrunken centroid fit.
#'
#' @references T. Hastie, R. Tibshirani, Balasubramanian Narasimhan and Gil Chu (2014).
#' pamr: Pam: prediction analysis for microarrays. R package version 1.55.
#' https://CRAN.R-project.org/package=pamr
#' @param pam.intcv.model a PAM classifier built with \code{pam.intcv}.
#' @param pred.obj dataset to have its sample group predicted.
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
#' biological.effect <- estimate.biological.effect(uhdata = uhdata.pl)
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#' biological.effect.nc <- biological.effect[!rownames(biological.effect) %in% ctrl.genes, ]
#' group.id <- substr(colnames(biological.effect.nc), 7, 7)
#'
#' biological.effect.train.ind <- colnames(biological.effect.nc)[c(sample(which(group.id == "E"), size = 64),
#'                                           sample(which(group.id == "V"), size = 64))]
#' biological.effect.test.ind <- colnames(biological.effect.nc)[!colnames(biological.effect.nc) %in% biological.effect.train.ind]
#'
#' biological.effect.nc.tr <- biological.effect.nc[, biological.effect.train.ind]
#' biological.effect.nc.te <- biological.effect.nc[, biological.effect.test.ind]
#'
#' pam.int <- pam.intcv(X = biological.effect.nc.tr,
#'                      y = substr(colnames(biological.effect.nc.tr), 7, 7),
#'                      kfold = 5, seed = 1)
#'
#' pam.pred <- pam.predict(pam.intcv.model = pam.int,
#'                         pred.obj = biological.effect.nc.te,
#'                         pred.obj.group.id = substr(colnames(biological.effect.nc.te), 7, 7))
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
