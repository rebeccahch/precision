#' Quantile nomalization
#'
#' Normalizes training dataset with quantile normalization, stores the quantiles from the training dataset as the references to frozen quantile normalize test dataset.
#' @references Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P. (2003) A Comparison of Normalization Methods for High Density Oligonucleotide Array Data Based on Bias and Variance. Bioinformatics 19(2) ,pp 185-193. http://bmbolstad.com/misc/normalize/normalize.html
#' @param train training data to be quantile normalized, rows as probes, columns as samples.
#' @param test test data to be frozen quantile normalized, rows as probes with equal number of rows as the training set, columns as samples.
#' @return a list of two datasets, the normalized training set and the frozen normalied test set
#' @export
#' @keywords preprocess
#' @import preprocessCore
#' @examples
#' set.seed(101)
#' group.id <- substr(colnames(non.r.data.pl), 7, 7)
#' train.ind <- colnames(non.r.data.pl)[c(sample(which(group.id == "E"), size = 64),
#'                                sample(which(group.id == "V"), size = 64))]
#' train.dat <- non.r.data.pl[, train.ind]
#' test.dat <- non.r.data.pl[, !colnames(non.r.data.pl) %in%train.ind]
#' data.qn <- quant.norm(train = train.dat)
#' str(data.qn)
#' data.qn <- quant.norm(train = train.dat, test = test.dat)
#' str(data.qn)


"quant.norm" <- function(train, test = NULL){
  stopifnot(nrow(train) == nrow(test))

  # quantile normalization training
  train.qn <- preprocessCore::normalize.quantiles(train, copy = TRUE)
  dimnames(train.qn) <- dimnames(train)

  # check if QN overrides train
  if(length(setdiff(train, train.qn)) == 0) stop("override train")

  if(is.null(test)){
    test.fqn <- NULL
  } else{
    #   record ref.dis
    ref.dis <- as.numeric(sort(train.qn[, 1]))

    ## frozen quantile normalize test
    test.fqn <- apply(test, 2, function(x){ord <- rank(x); ref.dis[ord]})
    dimnames(test.fqn) <- dimnames(test)
  }

  return(list("train.qn" = train.qn,
              "test.fqn" = test.fqn))
}
