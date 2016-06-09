#' Prediction with least absolute shrinkage and selection operator classifier
#'
#' Predicts from a least absolute shrinkage and selection operator fit.
#'
#' @param lasso.intcv.model a LASSO classifier built with lasso.intcv().
#' @param pred.obj expression dataset to have its sample group predicted, rows as probes, columns as samples; should have equal number of probes as the data trained.
#' @param pred.obj.group.id sample group corresponding to the dataset to be predicted; should have equal length as the number of samples as pred.obj.
#' @return predicted object, predicted error and predicted features
#' @export
#' @import glmnet
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
#' lasso.int <- lasso.intcv(X = smp.eff.nc.tr,
#'                          y = substr(colnames(smp.eff.nc.tr), 7, 7),
#'                          kfold = 5, seed = 1, alp = 1)
#'
#' lasso.pred <- lasso.predict(lasso.intcv.model = lasso.int,
#'                             pred.obj = smp.eff.nc.te,
#'                             pred.obj.group.id = substr(colnames(smp.eff.nc.te), 7, 7))
#' lasso.int$mc
#' lasso.pred$mc

"lasso.predict" <- function(lasso.intcv.model, pred.obj, pred.obj.group.id){
  pred <- predict(lasso.intcv.model$model, newx = t(pred.obj),
                  s = lasso.intcv.model$model$lambda.1se,
                  type = "class")

  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
  feature <- predict(lasso.intcv.model$model, newx = t(pred.obj),
                     s = lasso.intcv.model$model$lambda.1se)
  return(list(pred=pred, mc=mc, feature=feature))
}
