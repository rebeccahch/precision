#' Prediction with nearest shrunken centroid classifier
#'
#' Predict from a nearest shrunken centroid fit.
#'
#' @references T. Hastie, R. Tibshirani, Balasubramanian Narasimhan and Gil Chu (2014).
#' pamr: Pam: prediction analysis for microarrays. R package version 1.55.
#' https://CRAN.R-project.org/package=pamr
#' @param pam.intcv.model a PAM classifier built with \code{pam.intcv()}.
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
#' sample.effect <- estimate.sample.effect(uhdata = uhdata.pl)
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#' sample.effect.nc <- sample.effect[!rownames(sample.effect) %in% ctrl.genes, ]
#' group.id <- substr(colnames(sample.effect.nc), 7, 7)
#'
#' sample.effect.train.ind <- colnames(sample.effect.nc)[c(sample(which(group.id == "E"), size = 64),
#'                                           sample(which(group.id == "V"), size = 64))]
#' sample.effect.test.ind <- colnames(sample.effect.nc)[!colnames(sample.effect.nc) %in% sample.effect.train.ind]
#'
#' sample.effect.nc.tr <- sample.effect.nc[, sample.effect.train.ind]
#' sample.effect.nc.te <- sample.effect.nc[, sample.effect.test.ind]
#'
#' pam.int <- pam.intcv(X = sample.effect.nc.tr,
#'                      y = substr(colnames(sample.effect.nc.tr), 7, 7),
#'                      kfold = 5, seed = 1)
#'
#' pam.pred <- pam.predict(pam.intcv.model = pam.int,
#'                         pred.obj = sample.effect.nc.te,
#'                         pred.obj.group.id = substr(colnames(sample.effect.nc.te), 7, 7))
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
