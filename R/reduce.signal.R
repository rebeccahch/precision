#' Biological signal reduction
#'
#' Reduce biological signal by decreasing the mean group difference between sample groups.
#'
#' @param biological.effect the estimated biological effect dataset. The dataset must have rows as probes and columns as samples.
#' It can only take in probe-level dataset with a fixed number of probes per unique probe-set.
#' @param group.id a vector of sample-group labels for each sample of the estimated biological effect dataset.
#' @param group.id.level a vector of sample-group label level. It must have two and only two elements and the first element is the reference.
#' By default, \code{group.id.level = c("E", "V")}. That is in our study, we compare endometrial tumor samples to
#' ovarian tumor samples, with endometrial as our reference.
#' @param reduce.multiplier a multiplier specified to reduce between-sample-group signal by. By default, \code{reduce.multiplier = 1/2}.
#' @param pbset.id a vector of unique probe-set names. If it is not specified, it is the unique probe names of the dataset,
#' extracting from the row names.
#' @return estimated biological effect data, with reduced biological signal
#' @keywords data.setup
#' @export
#' @examples
#' biological.effect <- estimate.biological.effect(uhdata = uhdata.pl)
#' handling.effect <- estimate.handling.effect(uhdata = uhdata.pl,
#'                              nuhdata = nuhdata.pl)
#'
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#'
#' biological.effect.nc <- biological.effect[!rownames(biological.effect) %in% ctrl.genes, ]
#' handling.effect.nc <- handling.effect[!rownames(handling.effect) %in% ctrl.genes, ]
#' group.id <- substr(colnames(biological.effect.nc), 7, 7)
#'
#' redhalf.biological.effect.nc <- reduce.signal(biological.effect = biological.effect.nc,
#'                                     group.id = group.id,
#'                                     group.id.level = c("E", "V"),
#'                                     reduce.multiplier = 1/2)

"reduce.signal" <- function(biological.effect,
                                 group.id,
                                 group.id.level = c("E", "V"),
                                 reduce.multiplier = 1/2,
                                 pbset.id = NULL){
  stopifnot(nrow(biological.effect) != length(unique(rownames(biological.effect)))) # probe level
  stopifnot(length(unique(table(rownames(biological.effect)))) == 1)
  if(is.null(pbset.id)) pbset.id <- unique(rownames(biological.effect))

  n.p.u <- unique(table(rownames(biological.effect)))
  biological.effect.psl <- med.sum.pbset(biological.effect, num.per.unipbset = n.p.u)
  s.e.limma <- limma.pbset(data = biological.effect.psl,
                           group.id = group.id,
                           group.id.level = group.id.level,
                           pbset.id = pbset.id)
  de.ind <- s.e.limma$P.Value < 0.01

  sample.g1 <- rowMeans(biological.effect[rownames(biological.effect) %in% pbset.id[de.ind],
                                group.id == group.id.level[1]])
  sample.g2 <- rowMeans(biological.effect[rownames(biological.effect) %in% pbset.id[de.ind],
                                group.id == group.id.level[2]])

  half.signal <- (sample.g1 - sample.g2)*reduce.multiplier

  reduced.biological.effect.de <- cbind(biological.effect[rownames(biological.effect) %in% pbset.id[de.ind],
                                  group.id == group.id.level[1]] - half.signal,
                             biological.effect[rownames(biological.effect) %in% pbset.id[de.ind],
                                  group.id == group.id.level[2]])

  # combine and colnames, rownames back to original order
  temp <- biological.effect

  temp[rownames(temp) %in% pbset.id[de.ind], ] <- reduced.biological.effect.de[, colnames(biological.effect)]
  redhalf.biological.effect.pl.p10 <- temp; rm(temp)

  return(redhalf.biological.effect.pl.p10)
}
