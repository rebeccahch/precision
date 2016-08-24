#' Biological signal reduction
#'
#' Reduce biological signal by decreasing the mean group difference between sample groups.
#'
#' @param sample.effect the estimated sample effect dataset. The dataset must have rows as probes and columns as samples.
#' It can only take in probe-level dataset with a fixed number of probes per unique probe-set.
#' @param group.id a vector of sample-group labels for each sample of the estimated sample effect dataset.
#' @param group.id.level a vector of sample-group label level. It must have two and only two elements and the first element is the reference.
#' By default, \code{group.id.level = c("E", "V")}. That is in our study, we compare endometrial tumor samples to
#' ovarian tumor samples, with endometrial as our reference.
#' @param reduce.multiplier a multiplier specified to reduce between-sample-group signal by. By default, \code{reduce.multiplier = 1/2}.
#' @param pbset.id a vector of unique probe-set names. If it is not specified, it is the unique probe names of the dataset,
#' extracting from the row names.
#' @return estimated sample effect data, with reduced biological signal
#' @keywords data.setup
#' @export
#' @examples
#' sample.effect <- estimate.sample.effect(uhdata = uhdata.pl)
#' array.effect <- estimate.array.effect(uhdata = uhdata.pl,
#'                              nuhdata = nuhdata.pl)
#'
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#'
#' sample.effect.nc <- sample.effect[!rownames(sample.effect) %in% ctrl.genes, ]
#' array.effect.nc <- array.effect[!rownames(array.effect) %in% ctrl.genes, ]
#' group.id <- substr(colnames(sample.effect.nc), 7, 7)
#'
#' redhalf.sample.effect.nc <- reduce.signal(sample.effect = sample.effect.nc,
#'                                     group.id = group.id,
#'                                     group.id.level = c("E", "V"),
#'                                     reduce.multiplier = 1/2)

"reduce.signal" <- function(sample.effect,
                                 group.id,
                                 group.id.level = c("E", "V"),
                                 reduce.multiplier = 1/2,
                                 pbset.id = NULL){
  stopifnot(nrow(sample.effect) != length(unique(rownames(sample.effect)))) # probe level
  stopifnot(length(unique(table(rownames(sample.effect)))) == 1)
  if(is.null(pbset.id)) pbset.id <- unique(rownames(sample.effect))

  n.p.u <- unique(table(rownames(sample.effect)))
  sample.effect.psl <- med.sum.pbset(sample.effect, num.per.unipbset = n.p.u)
  s.e.limma <- limma.pbset(data = sample.effect.psl,
                           group.id = group.id,
                           group.id.level = group.id.level,
                           pbset.id = pbset.id)
  de.ind <- s.e.limma$P.Value < 0.01

  sample.g1 <- rowMeans(sample.effect[rownames(sample.effect) %in% pbset.id[de.ind],
                                group.id == group.id.level[1]])
  sample.g2 <- rowMeans(sample.effect[rownames(sample.effect) %in% pbset.id[de.ind],
                                group.id == group.id.level[2]])

  half.signal <- (sample.g1 - sample.g2)*reduce.multiplier

  reduced.sample.effect.de <- cbind(sample.effect[rownames(sample.effect) %in% pbset.id[de.ind],
                                  group.id == group.id.level[1]] - half.signal,
                             sample.effect[rownames(sample.effect) %in% pbset.id[de.ind],
                                  group.id == group.id.level[2]])

  # combine and colnames, rownames back to original order
  temp <- sample.effect

  temp[rownames(temp) %in% pbset.id[de.ind], ] <- reduced.sample.effect.de[, colnames(sample.effect)]
  redhalf.sample.effect.pl.p10 <- temp; rm(temp)

  return(redhalf.sample.effect.pl.p10)
}
