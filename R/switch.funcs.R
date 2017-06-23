"switch.norm.funcs" <- function(norm.list = c("NN", "QN"),
                                 norm.funcs = NULL){
  temp <- tolower(norm.list) %in% c("nn", "mn", "qn", "vsn")
  stopifnot(length(unique(norm.list[!temp])) == length(norm.funcs))

  new.norm.funcs <- tolower(norm.list)
  new.norm.funcs <- gsub("nn", "nn.norm",
                         gsub("mn", "med.norm",
                              gsub("qn", "quant.norm",
                                   gsub("vsn", "vs.norm",
                                        tolower(new.norm.funcs)))))
  new.norm.funcs <- c(new.norm.funcs[temp], norm.funcs)
  return(new.norm.funcs)
}

"switch.norm.funcs.flex" <- function(norm.list = c("NN", "QN"),
                                norm.funcs = NULL){

  temp <- tolower(norm.list) %in% c("nn", "mn", "qn", "vsn")
  stopifnot(length(unique(norm.list[!temp])) == length(norm.funcs))

  new.norm.funcs <- tolower(norm.list)

  switch_x <- sapply(new.norm.funcs, function(x) switch(x,
                                                        "nn" = "nn.norm",
                                                        "mn" = "med.norm",
                                                        "qn" = "quant.norm",
                                                        "vsn" = "vs.norm"))

  switch_x <- as.character(switch_x)
  switch_x[which(switch_x == "NULL")] <-
    paste0(tolower(norm.list)[which(switch_x == "NULL")], ".norm")
  new.norm.funcs <- switch_x

  stopifnot(new.norm.funcs %in% c("nn.norm", "med.norm", "quant.norm", "vs.norm",
                                  norm.funcs))

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
