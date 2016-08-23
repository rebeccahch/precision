#' Median normalization
#'
#' Normalize a training dataset so that each array shares
#' a same median and store the median from the training dataset
#' as the reference to frozen median normalize a test dataset.
#' Also two other options are available: to only normalize
#' a training dataset but not frozen normalize a test dataset, or vise versa.
#'
#' @param train the training dataset to be median normalized.
#' The dataset must have rows as probes and columns as samples.
#' This can be left unspecified if \code{ref.dis} is
#' suppied for frozen normalize test set.
#' @param test the test dataset to be frozen median normalized.
#' The dataset must have rows as probes and columns as samples.
#' The number of rows must equal to the number of rows in the training set.
#' By default, the test set is not specified (\code{test = NULL}) and
#' no frozen normalization will be performed.
#' @param ref.dis the reference distribution for frozen median normalize test set
#' against previously normalized training set.
#' This is required when \code{train} is not supplied.
#' By default, \code{ref.dis = NULL}.
#' @return a list of two datasets and one reference distribution:
#' \item{train.mn}{the normalized training set, if training set is specified}
#' \item{test.fmn}{the frozen normalized test set, if test set is specified}
#' \item{ref.dis}{the reference distribution}
#' @importFrom stats median
#' @export
#' @keywords preprocess
#' @examples
#' set.seed(101)
#' group.id <- substr(colnames(nuhdata.pl), 7, 7)
#' train.ind <- colnames(nuhdata.pl)[c(sample(which(group.id == "E"), size = 64),
#'                                sample(which(group.id == "V"), size = 64))]
#' train.dat <- nuhdata.pl[, train.ind]
#' test.dat <- nuhdata.pl[, !colnames(nuhdata.pl) %in% train.ind]
#'
#' # normalize only training set
#' data.mn <- med.norm(train = train.dat)
#' str(data.mn)
#'
#' # normalize training set and frozen normalize test set
#' data.mn <- med.norm(train = train.dat, test = test.dat)
#' str(data.mn)
#'
#' # frozen normalize test set with reference distribution
#' ref <- med.norm(train = train.dat)$ref.dis
#' data.mn <- med.norm(test = test.dat, ref.dis = ref)
#' str(data.mn)
#'

"med.norm" <- function(train = NULL, test = NULL, ref.dis = NULL){
  if(!is.null(train) && !is.null(test)) stopifnot(nrow(train) == nrow(test))
  if(is.null(train)) stopifnot(!is.null(ref.dis))

  if(!is.null(train)){
    # median normalization training
    ref.dis <- median(train)
    temp <- apply(train, 2, median) - ref.dis
    shifts.train <- matrix(rep(temp, each = nrow(train)), ncol = ncol(train))
    train.mn <- train - shifts.train
  } else train.mn <- NULL

  if(is.null(test)) {
    test.fmn <- NULL
  } else{
    temp <- apply(test, 2, median) - ref.dis
    shifts.test <- matrix(rep(temp, each = nrow(test)), ncol = ncol(test))
    test.fmn <- test - shifts.test
  }

  return(list("train.mn" = train.mn,
              "test.fmn" = test.fmn,
              "ref.dis" = ref.dis))
}
