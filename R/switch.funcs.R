"switch.norm.funcs" <- function(norm.list = c("NN", "QN"),
                                norm.funcs = NULL){
  temp <- tolower(norm.list) %in% c("nn", "mn", "qn", "vsn")
  stopifnot(length(norm.list[!temp]) == length(norm.funcs))

  new.norm.funcs <- tolower(norm.list)
  new.norm.funcs <- gsub("mn", "med.norm",
                         gsub("qn", "quant.norm",
                              gsub("vsn", "vs.norm", tolower(new.norm.funcs))))
  new.norm.funcs <- c(new.norm.funcs[temp], norm.funcs)
  return(new.norm.funcs)
}


"switch.classifier.funcs" <- function(class.list = c("PAM", "LASSO"),
                                      class.funcs = NULL,
                                      pred.funcs = NULL){
  temp <- tolower(class.list) %in% c("pam", "lasso")
  stopifnot(length(class.list[!temp]) == length(class.funcs))
  stopifnot(length(class.list[!temp]) == length(pred.funcs))

  new.class.funcs <- tolower(class.list)
  new.class.funcs <- gsub("pam", "pam.intcv",
                         gsub("lasso", "lasso.intcv", tolower(new.class.funcs)))
  new.class.funcs <- c(new.class.funcs[temp], class.funcs)

  new.pred.funcs <- tolower(class.list)
  new.pred.funcs <- gsub("pam", "pam.predict",
                          gsub("lasso", "lasso.predict", tolower(new.pred.funcs)))
  new.pred.funcs <- c(new.pred.funcs[temp], pred.funcs)

  return(list(build.funcs = new.class.funcs,
              pred.funcs = new.pred.funcs))
}
