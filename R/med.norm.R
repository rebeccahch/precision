#' Median nomalization
#'
#' Normalizes training dataset so that each array shares a same median, stores the median from the training dataset as the reference to frozen median normalize test dataset.
#'
#' @param train training data to be median normalized, rows as probes, columns as samples.
#' @param test test data to be frozen median normalized, rows as probes with equal number of rows as the training set, columns as samples.
#' @return a list of two datasets, the normalized training set and the frozen normalied test set
#' @importFrom stats median
#' @export
#' @keywords preprocess
#' @examples
#' set.seed(101)
#' group.id <- substr(colnames(non.r.data.pl), 7, 7)
#' train.ind <- colnames(non.r.data.pl)[c(sample(which(group.id == "E"), size = 64),
#'                                sample(which(group.id == "V"), size = 64))]
#' train.dat <- non.r.data.pl[, train.ind]
#' test.dat <- non.r.data.pl[, !colnames(non.r.data.pl) %in%train.ind]
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
