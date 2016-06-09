"create.storage" <- function(design.list = c("CC+", "CC-", "PC+", "PC-"),
                             class.list = c("PAM", "LASSO"),
                             norm.list = c("NN", "QN"),
                             validating.sets = c("test", "gold.standard")){
  n.design <- length(design.list)
  n.norm <- length(norm.list)
  n.class <- length(class.list)
  n.valid <- length(validating.sets)

  storage <- matrix(rep(list(rep(list(NA), n.norm*n.class)), n.design*n.valid),
                    ncol = n.design, nrow = n.valid)
  colnames(storage) <- design.list
  rownames(storage) <- validating.sets
  for(i in 1:nrow(storage)){
    for(j in 1:ncol(storage)){
      names(storage[i, j][[1]]) <- paste(rep(class.list, each = n.norm), norm.list, sep = ".")
    }
  }

  return(storage)
}
