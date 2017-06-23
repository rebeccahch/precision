#' Nearest shrunken centroid through internal cross validation
#'
#' Build a PAM classifier using internal cross validation to choose
#' the tuning parameter, with 5-fold cross validation as the default.
#'
#' @references T. Hastie, R. Tibshirani, Balasubramanian Narasimhan and Gil Chu (2014).
#' pamr: Pam: prediction analysis for microarrays. R package version 1.55.
#' https://CRAN.R-project.org/package=pamr
#' @param X dataset to be trained.
#' This dataset must have rows as probes and columns as samples.
#' @param y a vector of sample group of each sample for the dataset to be trained.
#' It must have an equal length to the number of samples in \code{X}.
#' @param vt.k custom-specified threshold list.
#' By default, \code{vt.k = NULL} and
#' 30 values will be predetermined by the pamr package.
#' @param n.k number of threshold values desired. By default, \code{n.k = 30}.
#' @param kfold number of folds. By default, \code{kfold = 5}.
#' @param folds pre-specifies samples to each fold.
#' By default, \code{folds = NULL} for no pre-specification.
#' @param seed an integer used to initialize a pseudorandom number generator.
#' @return a list of 4 elements:
#' \item{mc}{an internal misclassification error rate}
#' \item{time}{processing time of performing internal validation with PAM}
#' \item{model}{a PAM classifier, resulted from \code{pamr.train}}
#' \item{cfs}{estimated coefficients for the final classifier}
#' @import pamr
#' @export
#' @keywords classification
#' @examples
#' set.seed(101)
#' biological.effect <- estimate.biological.effect(uhdata = uhdata.pl)
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#' biological.effect.nc <- biological.effect[!rownames(biological.effect) %in% ctrl.genes, ]
#' group.id <- substr(colnames(biological.effect.nc), 7, 7)
#'
#' biological.effect.train.ind <- colnames(biological.effect.nc)[c(sample(which(group.id == "E"), size = 64),
#'                                          sample(which(group.id == "V"), size = 64))]
#' biological.effect.nc.tr <- biological.effect.nc[, biological.effect.train.ind]
#'
#' pam.int <- pam.intcv(X = biological.effect.nc.tr,
#'                      y = substr(colnames(biological.effect.nc.tr), 7, 7),
#'                      kfold = 5, seed = 1)
#'

"pam.intcv" <- function(X, y, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL, seed){

  ptm <- proc.time()
  set.seed(seed)
  data.pam  <- list(x = X, y = factor(y), geneids = rownames(X), genenames = rownames(X))
  fit.pam	<- pamr::pamr.train(data.pam, threshold=vt.k, n.threshold = n.k)
  fit.cv <-  new.pamr.cv(fit = fit.pam, data = data.pam, nfold = kfold)
  best.threshold <- fit.cv$threshold[max(which(fit.cv$error == min(fit.cv$error)))]

  mc <- fit.cv$error[which.min(fit.cv$error)]

  model <- pamr::pamr.train(data.pam, threshold = best.threshold, n.threshold = n.k)

  ## if nonzero == 0 (no feature selected)
  coefs <- trycatch.func(pamr::pamr.listgenes(model, data.pam, threshold = best.threshold))

  time <- proc.time() - ptm
  return(list(mc = mc, time = time, model = model, cfs = coefs))
}



"new.pamr.cv" <- function (fit, data, nfold = 5, ...){
  x <- data$x[fit$gene.subset, fit$sample.subset]
  if (is.null(fit$newy)) {
    y <- factor(data$y[fit$sample.subset])
  }
  else {
    y <- factor(data$newy[fit$sample.subset])
  }
  this.call <- match.call()
  nsccv2 <- get("nsccv", envir = asNamespace("pamr"))
  balanced.folds <- get("balanced.folds", envir = asNamespace("pamr"))
  folds = balanced.folds(y, nfolds = nfold)
  junk <- nsccv2(x, y, object = fit, folds = folds,
                 survival.time = data$survival.time,
                 censoring.status = data$censoring.status,
                 ngroup.survival = fit$ngroup.survival,
                 problem.type = fit$problem.type,
                 ...) # changed here
  junk$call <- this.call
  junk$newy <- fit$newy
  junk$sample.subset <- fit$sample.subset
  return(junk)
}
