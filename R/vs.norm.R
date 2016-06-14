#' Variance stabilizing normalization
#'
#' Normalize training dataset with vsn and
#' store the fitted vsn model from the training dataset as the reference to frozen variance stabilizing normalize test dataset.
#'
#' @references Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, Annemarie Poustka and Martin Vingron.
#' Variance Stabilization Applied to Microarray Data Calibration and to the Quantification of Differential Expression.
#' Bioinformatics 18, S96-S104 (2002).
#' @param train training dataset to be variance stabilizing normalized. The dataset must have rows as probes and columns as samples.
#' @param test test dataset to be frozen variance stabilizing normalized. The dataset must have rows as probes and columns as samples.
#' The number of rows must equal to the number of rows in the training set.
#' By default, the test set is not specified (\code{test = NULL}) and no frozen normalization will be performed.
#' @return a list of two datasets:
#' \item{train.mn}{the normalized training set}
#' \item{test.fmn}{the frozen normalized test set, if test set is specified}
#' @export
#' @keywords preprocess
#' @import vsn
#' @examples
#' \dontrun{
#' set.seed(101)
#' group.id <- substr(colnames(non.r.data.pl), 7, 7)
#' train.ind <- colnames(non.r.data.pl)[c(sample(which(group.id == "E"), size = 64),
#'                                sample(which(group.id == "V"), size = 64))]
#' train.dat <- non.r.data.pl[, train.ind]
#' test.dat <- non.r.data.pl[, !colnames(non.r.data.pl) %in%train.ind]
#' data.vsn <- vs.norm(train = train.dat)
#' str(data.vsn)
#' data.vsn <- vs.norm(train = train.dat, test = test.dat)
#' str(data.vsn)
#' }

"vs.norm" <- function(train, test = NULL){
  stopifnot(nrow(train) == nrow(test))

  # train -> train.vsn
  train0 <- 2^train
  train.vsn0 <- vsn::vsn2(train0)
  train.vsn <- log2(exp(as.matrix(train.vsn0)))

  if(is.null(test)) {
    test.fvsn <- NULL
  } else {
    # test -> test.fvsn w.r.t. train.vsn
    test0 <- 2^test
    test.fvsn0 <- vsn::vsn2(test0, train.vsn0)
    test.fvsn <- log2(exp(as.matrix(test.fvsn0)))
  }

  return(list("train.vsn" = train.vsn,
              "test.fvsn" = test.fvsn))
}
