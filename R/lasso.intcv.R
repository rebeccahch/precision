#' Least absolute shrinkage and selection operator through internal cross validation
#'
#' Builds a LASSO classifier using internal cross validation, with a 5-fold cross validation as default.
#'
#' @param X expression dataset to be trained, rows as probes, columns as samples.
#' @param y sample group corresponding to the dataset to be trained; should have the equal length as the number of samples as X.
#' @param kfold number of folds; default is 5.
#' @param seed specifies seed for random assignment using set.seed().
#' @param alp alpha, the penalty type from 0 to 1; default alp = 1 for LASSO; alp = 0 for a ridge classifier.
#' @return a LASSO classifier
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
#'                                          sample(which(group.id == "V"), size = 64))]
#' smp.eff.nc.tr <- smp.eff.nc[, smp.eff.train.ind]
#'
#' lasso.int <- lasso.intcv(X = smp.eff.nc.tr,
#'                          y = substr(colnames(smp.eff.nc.tr), 7, 7),
#'                          kfold = 5, seed = 1, alp = 1)

"lasso.intcv" <- function(kfold = 5, X, y, seed, alp = 1){
  ptm <- proc.time()
  set.seed(seed)

  cv.fit <- glmnet::cv.glmnet(x = data.matrix(t(X)), y = factor(y),
                      family = "binomial", type.measure = "class", alpha = alp, nfold = kfold)
  mc <- cv.fit$cvm[which(cv.fit$lambda == cv.fit$lambda.1se)]
  #best.lambda <- cv.fit$lambda.1se # can be extracted from cv.fit
  coefs <- trycatch.func(coef(cv.fit, s = "lambda.1se"))
  time <- proc.time() - ptm
  return(list(mc=mc, time=time, model=cv.fit, cfs=coefs))
}
