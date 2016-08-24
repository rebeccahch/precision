#' Variance stabilizing normalization
#'
#' Normalize training dataset with vsn and
#' store the fitted vsn model from the training dataset as the reference
#' to frozen variance stabilizing normalize test dataset.
#' Also two other options are available: to only normalize
#' a training dataset but not frozen normalize a test dataset, or vise versa.
#'
#' @references Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, Annemarie Poustka and Martin Vingron.
#' Variance Stabilization Applied to Microarray Data Calibration and to the Quantification of Differential Expression.
#' Bioinformatics 18, S96-S104 (2002).
#' @param train training dataset to be variance stabilizing normalized.
#' The dataset must have rows as probes and columns as samples.
#' This can be left unspecified if \code{ref.dis} is
#' suppied for frozen normalize test set.
#' @param test test dataset to be frozen variance stabilizing normalized.
#' The dataset must have rows as probes and columns as samples.
#' The number of rows must equal to the number of rows in the training set.
#' By default, the test set is not specified (\code{test = NULL}) and
#' no frozen normalization will be performed.
#' @param ref.dis reference distribution for frozen variance stabilizing normalize test set
#' against previously normalized training set. This is required when \code{train} is not supplied.
#' By default, \code{ref.dis = NULL}.
#' @return a list of two datasets and one reference distribution:
#' \item{train.mn}{the normalized training set}
#' \item{test.fmn}{the frozen normalized test set, if test set is specified}
#' \item{ref.dis}{the reference distribution}
#' @export
#' @keywords preprocess
#' @import vsn
#' @examples
#' \dontrun{
#' set.seed(101)
#' group.id <- substr(colnames(nuhdata.pl), 7, 7)
#' train.ind <- colnames(nuhdata.pl)[c(sample(which(group.id == "E"), size = 64),
#'                                sample(which(group.id == "V"), size = 64))]
#' train.dat <- nuhdata.pl[, train.ind]
#' test.dat <- nuhdata.pl[, !colnames(nuhdata.pl) %in% train.ind]
#'
#' # normalize only training set
#' data.vsn <- vs.norm(train = train.dat)
#' str(data.vsn)
#'
#' # normalize training set and frozen normalize test set
#' data.vsn <- vs.norm(train = train.dat, test = test.dat)
#' str(data.vsn)
#'
#' # frozen normalize test set with reference distribution
#' ref <- vs.norm(train = train.dat)$ref.dis
#' data.vsn <- vs.norm(test = test.dat, ref.dis = ref)
#' str(data.vsn)
#' }
#'

"vs.norm" <- function(train = NULL, test = NULL, ref.dis = NULL){
  if(!is.null(train) && !is.null(test)) stopifnot(nrow(train) == nrow(test))
  if(is.null(train)) stopifnot(!is.null(ref.dis))

  if(!is.null(train)){
    # vsn training
    train0 <- 2^train
    ref.dis <- vsn::vsn2(train0)
    train.vsn <- log2(exp(as.matrix(ref.dis)))
  } else train.vsn <- NULL

  if(is.null(test)) {
    test.fvsn <- NULL
  } else {
    test0 <- 2^test
    test.fvsn0 <- vsn::vsn2(test0, ref.dis)
    test.fvsn <- log2(exp(as.matrix(test.fvsn0)))
  }

  return(list("train.vsn" = train.vsn,
              "test.fvsn" = test.fvsn,
              "ref.dis" = ref.dis))
}
