#' Nearest shrunken centroid through internal cross validation
#'
#' Builds a PAM classifier using internal cross validation, with 5-fold cross validation as the default.
#'
#' @param X expression dataset to be trained, rows as probes, columns as samples.
#' @param y sample group corresponding to the data to be trained; must have the equal length as the number of samples as X.
#' @param vt.k custom-specified threshold list; default is NULL predetermined by the PAM package.
#' @param n.k number of threshold values desired; default is 30.
#' @param kfold number of folds for cross validation; default is 5
#' @param folds prespecifies samples to folds; default is NULL for no prespecification.
#' @param seed specifies seed for random assignment using set.seed().
#' @return a PAM classifier
#' @import pamr
#' @export
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
#' pam.int <- pam.intcv(X = smp.eff.nc.tr,
#'                      y = substr(colnames(smp.eff.nc.tr), 7, 7),
#'                      kfold = 5, seed = 1)

"pam.intcv" <- function(X, y, vt.k=NULL, n.k=30, kfold = 5, folds=NULL, seed){

  ptm <- proc.time()
  set.seed(seed)
  data.pam  <- list(x=X, y=factor(y), geneids = rownames(X), genenames = rownames(X))
  fit.pam	<- pamr::pamr.train(data.pam, threshold=vt.k, n.threshold=n.k)
  fit.cv <-  new.pamr.cv(fit=fit.pam, data=data.pam, nfold = kfold)
  best.threshold <- fit.cv$threshold[max(which(fit.cv$error==min(fit.cv$error)))]

  mc <- fit.cv$error[which.min(fit.cv$error)]

  model <- pamr::pamr.train(data.pam, threshold=best.threshold, n.threshold=n.k)

  ## if nonzero == 0 (no feature selected)
  coefs <- trycatch.func(pamr::pamr.listgenes(model, data.pam, threshold = best.threshold))

  time <- proc.time() - ptm
  return(list(mc=mc, time=time, model = model, cfs = coefs))
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
  junk <- nsccv2(x, y, object = fit, folds = folds, survival.time = data$survival.time, censoring.status = data$censoring.status,
                 ngroup.survival = fit$ngroup.survival, problem.type = fit$problem.type,
                 ...) # changed here
  junk$call <- this.call
  junk$newy <- fit$newy
  junk$sample.subset <- fit$sample.subset
  return(junk)
}
