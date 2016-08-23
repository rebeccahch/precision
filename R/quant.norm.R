#' Quantile normalization
#'
#' Normalize training dataset with quantile normalization and
#' store the quantiles from the training dataset as the references to frozen quantile normalize test dataset.
#'
#' @references Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P. (2003)
#' A Comparison of Normalization Methods for High Density Oligonucleotide Array Data Based on Bias and Variance.
#' Bioinformatics 19(2) , pp 185-193.
#' http://bmbolstad.com/misc/normalize/normalize.html
#' @param train training dataset to be quantile normalized.
#' The dataset must have rows as probes and columns as samples.
#' This can be left unspecified if \code{ref.dis} is
#' suppied for frozen normalize test set.
#' @param test test dataset to be frozen quantile normalized.
#' The dataset must have rows as probes and columns as samples.
#' The number of rows must equal to the number of rows in the training set.
#' By default, the test set is not specified (\code{test = NULL}) and no frozen normalization will be performed.
#' @param ref.dis reference distribution for frozen quantile normalize test set
#' against previously normalized training set. This is required when \code{train} is not supplied.
#' By default, \code{ref.dis = NULL}.
#' @return a list of two datasets and one reference distribution:
#' \item{train.mn}{the normalized training set}
#' \item{test.fmn}{the frozen normalized test set, if test set is specified}
#' \item{ref.dis}{the reference distribution}
#' @export
#' @keywords preprocess
#' @import preprocessCore
#' @examples
#' set.seed(101)
#' group.id <- substr(colnames(nuhdata.pl), 7, 7)
#' train.ind <- colnames(nuhdata.pl)[c(sample(which(group.id == "E"), size = 64),
#'                                sample(which(group.id == "V"), size = 64))]
#' train.dat <- nuhdata.pl[, train.ind]
#' test.dat <- nuhdata.pl[, !colnames(nuhdata.pl) %in% train.ind]
#'
#' # normalize only training set
#' data.qn <- quant.norm(train = train.dat)
#' str(data.qn)
#'
#' # normalize training set and frozen normalize test set
#' data.qn <- quant.norm(train = train.dat, test = test.dat)
#' str(data.qn)
#'
#' # frozen normalize test set with reference distribution
#' ref <- quant.norm(train = train.dat)$ref.dis
#' data.qn <- quant.norm(test = test.dat, ref.dis = ref)
#' str(data.qn)
#'

"quant.norm" <- function(train = NULL, test = NULL, ref.dis = NULL){
  if(!is.null(train) && !is.null(test)) stopifnot(nrow(train) == nrow(test))
  if(is.null(train)) stopifnot(!is.null(ref.dis))

  if(!is.null(train)){
    # quantile normalization training
    train.qn <- preprocessCore::normalize.quantiles(train, copy = TRUE)
    dimnames(train.qn) <- dimnames(train)

    ref.dis <- as.numeric(sort(train.qn[, 1]))
  } else train.qn <- NULL

  if(is.null(test)){
    test.fqn <- NULL
  } else{
    ## frozen quantile normalize test
    test.fqn <- apply(test, 2, function(x){ord <- rank(x); ref.dis[ord]})
    dimnames(test.fqn) <- dimnames(test)
  }

  return(list("train.qn" = train.qn,
              "test.fqn" = test.fqn,
              "ref.dis" = ref.dis))
}
