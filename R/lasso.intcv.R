#' Least absolute shrinkage and selection operator through internal cross validation
#'
#' Build a LASSO classifier using internal cross validation to choose the turning parameter, with a 5-fold cross validation as default.
#'
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Mod- els via Coordinate Descent,
#' http://www.stanford.edu/~hastie/Papers/glmnet.pdf
#' Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010
#' @param X dataset to be trained. This dataset must have rows as probes and columns as samples.
#' @param y a vector of sample group of each sample for the dataset to be trained.
#' It must have an equal length to the number of samples in \code{X}.
#' @param kfold number of folds. By default, \code{kfold = 5}.
#' @param seed an integer used to initialize a pseudorandom number generator.
#' @param alp alpha, the penalty type. It can be any numeric value from 0 to 1.
#' By default, \code{alp = 1} which is for LASSO. \code{alp = 0} is for ridge and
#' any value in between is for elastic net.
#' @return a list of 4 elements:
#' \item{mc}{an internal misclassification error rate}
#' \item{time}{the processing time of performing internal validation with LASSO}
#' \item{model}{a LASSO classifier, resulted from \code{cv.fit}}
#' \item{cfs}{estimated coefficients for the final classifier}
#' @export
#' @import glmnet
#' @keywords classification
#' @examples
#' set.seed(101)
#' biological.effect <- estimate.biological.effect(uhdata = uhdata.pl)
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#' biological.effect.nc <- biological.effect[!rownames(biological.effect)
#'   %in% ctrl.genes, ]
#' group.id <- substr(colnames(biological.effect.nc), 7, 7)
#'
#' biological.effect.train.ind <- colnames(biological.effect.nc)[c(sample(which(
#'   group.id == "E"), size = 64),
#'   sample(which(group.id == "V"), size = 64))]
#' biological.effect.nc.tr <- biological.effect.nc[, biological.effect.train.ind]
#'
#' lasso.int <- lasso.intcv(X = biological.effect.nc.tr,
#'                          y = substr(colnames(biological.effect.nc.tr), 7, 7),
#'                          kfold = 5, seed = 1, alp = 1)
#'

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
