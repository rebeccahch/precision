#' Median normalization
#'
#' Normalize training dataset so that each array shares
#' a same median and store the median from the training dataset
#' as the reference to frozen median normalize test dataset.
#'
#' @param train training dataset to be median normalized.
#' The dataset must have rows as probes and columns as samples.
#' @param test test dataset to be frozen median normalized.
#' The dataset must have rows as probes and columns as samples.
#' The number of rows must equal to the number of rows in the training set.
#' By default, the test set is not specified (\code{test = NULL}) and
#' no frozen normalization will be performed.
#' @return a list of two datasets:
#' \item{train.mn}{the normalized training set}
#' \item{test.fmn}{the frozen normalized test set, if test set is specified}
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
#' data.mn <- med.norm(train = train.dat)
#' str(data.mn)
#' data.mn <- med.norm(train = train.dat, test = test.dat)
#' str(data.mn)

"med.norm" <- function(train, test = NULL){
  stopifnot(nrow(train) == nrow(test))

  # median normalization training
  temp <- apply(train, 2, median) - median(train)
  shifts.train <- matrix(rep(temp, each = nrow(train)), ncol = ncol(train))
  train.mn <- train - shifts.train
  #(train - train.mn) == shifts.train

  if(is.null(test)) {
    test.fmn <- NULL
  } else{
      # train median
      train.med <- median(train)

      # frozen quantile normalization test
      temp <- apply(test, 2, median) - train.med
      shifts.test <- matrix(rep(temp, each = nrow(test)), ncol = ncol(test))
      test.fmn <- test - shifts.test
  }

  return(list("train.mn" = train.mn,
              "test.fmn" = test.fmn))
}
