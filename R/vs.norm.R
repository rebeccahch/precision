#' Variance stabalizing nomalization
#'
#' Normalizes training dataset with vsn, stores the fitted vsn model from the training dataset as the reference to frozen variance stabalizng normalize test dataset.
#' @references Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, Annemarie Poustka and Martin Vingron. Variance Stabilization Applied to Microarray Data Calibration and to the Quantification of Differential Expression. Bioinformatics 18, S96-S104 (2002).
#' @param train training data to be variance stabalizing normalized, rows as probes, columns as samples.
#' @param test test data to be frozen variance stabalizing normalized, rows as probes with equal number of rows as the training set, columns as samples.
#' @return a list of two datasets, the normalized training set and the frozen normalied test set
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
