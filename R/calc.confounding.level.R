#' Level of confounding calculation
#'
#' Calculates level of confounding between handling effects of a dataset and sample-group labels.
#'
#' @param data expression dataset, rows as probes, columns as samples.
#' @param group.id sample group ID for the dataset.
#' @param nbe.genes a vector of non-biological genes indicated as TRUE or 1, equal length as number of probes in the data.
#' @return the level of confounding and the most correlated principal component of the non-biological genes in the dataset with the sample-group labels
#' @keywords data.setup
#' @importFrom stats lm prcomp
#' @export
#' @examples
#' smp.eff <- estimate.smp.eff(r.data = r.data.pl)
#' ary.eff <- estimate.ary.eff(r.data = r.data.pl,
#'                              non.r.data = non.r.data.pl)
#' ctrl.genes <- unique(rownames(r.data.pl))[grep("NC", unique(rownames(r.data.pl)))]
#' smp.eff.nc <- smp.eff[!rownames(smp.eff) %in% ctrl.genes, ]
#' ary.eff.nc <- ary.eff[!rownames(ary.eff) %in% ctrl.genes, ]
#' group.id <- substr(colnames(smp.eff.nc), 7, 7)
#' smp.eff.train.ind <- colnames(smp.eff.nc)[c(sample(which(group.id == "E"), size = 64),
#' sample(which(group.id == "V"), size = 64))]
#' ary.eff.train.ind <- colnames(ary.eff.nc)[c(1:64, 129:192)]
#' # randomly created a vector of Boolean for nbe.genes
#' nbe.genes <- sample(c(TRUE, FALSE), size = nrow(smp.eff.nc), replace = TRUE)
#' calc.confounding.level(data = smp.eff.nc[, smp.eff.train.ind],
#'                        group.id = substr(smp.eff.train.ind, 7, 7),
#'                        nbe.genes = nbe.genes)

"calc.confounding.level" <- function(data, group.id, nbe.genes){
  stopifnot(is.logical(nbe.genes))

  data.nbe <- data[nbe.genes, ]
  pca.data <- prcomp(t(data.nbe), scale = TRUE)

  temp <- rep(0,5)
  for(k in 1:5){
    temp[k] <- summary(lm(pca.data$x[,k] ~ group.id))$adj.r.squared
  }
  locc <- max(temp)
  k_pc <- which.max(temp)
  rm(temp)

  return(list(locc = locc, k_pc = k_pc))
}
