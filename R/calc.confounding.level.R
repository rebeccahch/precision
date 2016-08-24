#' Level of confounding calculation
#'
#' Calculate the level of confounding between handling effects and sample group of interest
#' for a dataset. First, principal component is applied on the non-biological subset of the data.
#' The first five principal components are then used to build a simple linear regression model to predict the sample group.
#' the highest adjusted R-squared is returned as the level of confounding.
#'
#' @references Leek J., Scharpf R., Bravo H., et al. Tackling the widespread and critical impact of batch effects in high-throughput data. Nat Rev Genet 11:733-9, 2010.
#' @param data microarry dataset. It must have rows as probes and columns as samples.
#' @param group.id a vector of sample-group labels for each sample of the dataset.
#' @param nbe.genes a vector of non-biological genes used to filter the dataset.
#' Non-biological genes are indicated as \code{TRUE}, otherwise as \code{FALSE}.
#' The vector must have an equal length to the number of probes in the dataset.
#' @return a list of two elements:
#' \item{locc}{the level of confounding}
#' \item{k_pc}{the most correlated principal component of the
#' non-biological genes in the dataset with the sample group}
#' @keywords data.setup
#' @importFrom stats lm prcomp
#' @export
#' @examples
#' \dontrun{
#' sample.effect <- estimate.sample.effect(uhdata = uhdata.pl)
#' array.effect <- estimate.array.effect(uhdata = uhdata.pl,
#'                              nuhdata = nuhdata.pl)
#'
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#'
#' sample.effect.nc <- sample.effect[!rownames(sample.effect) %in% ctrl.genes, ]
#' array.effect.nc <- array.effect[!rownames(array.effect) %in% ctrl.genes, ]
#'
#' group.id <- substr(colnames(sample.effect.nc), 7, 7)
#'
#' sample.effect.train.ind <- colnames(sample.effect.nc)[c(sample(which(group.id == "E"), size = 64),
#' sample(which(group.id == "V"), size = 64))]
#' array.effect.train.ind <- colnames(array.effect.nc)[c(1:64, 129:192)]
#'
#' # randomly created a vector of Boolean for nbe.genes
#' nbe.genes <- sample(c(TRUE, FALSE), size = nrow(sample.effect.nc), replace = TRUE)
#'
#' calc.confounding.level(data = sample.effect.nc[, sample.effect.train.ind],
#'                        group.id = substr(sample.effect.train.ind, 7, 7),
#'                        nbe.genes = nbe.genes)
#' }
#'

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
