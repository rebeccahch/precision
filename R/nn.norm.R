"nn.norm" <- function (train = NULL, test = NULL, ref.dis = NULL){
  if (!is.null(train) && !is.null(test)) 
    stopifnot(nrow(train) == nrow(test))
  if (is.null(train)) 
    stopifnot(!is.null(ref.dis))
  if (!is.null(train)){
    train.nn <- train
    dimnames(train.nn) <- dimnames(train)
    ref.dis <- list() # no need, but can't be NULL for stopifnot
  }
  else train.nn <- NULL
  if (is.null(test)){
    test.fnn <- NULL
  }
  else {
    test.fnn <- test
    dimnames(test.fnn) <- dimnames(test)
  }
  return(list(train.nn = train.nn, test.fnn = test.fnn, ref.dis = ref.dis))
}