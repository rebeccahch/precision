#' Biological signal reduction
#'
#' Reduce biological signal by decreasing the mean group difference between sample groups.
#'
#' @param smp.eff the estimated sample effect dataset. The dataset must have rows as probes and columns as samples.
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
#' smp.eff <- estimate.smp.eff(r.data = r.data.pl)
#' ary.eff <- estimate.ary.eff(r.data = r.data.pl,
#'                              non.r.data = non.r.data.pl)
#'
#' ctrl.genes <- unique(rownames(r.data.pl))[grep("NC", unique(rownames(r.data.pl)))]
#'
#' smp.eff.nc <- smp.eff[!rownames(smp.eff) %in% ctrl.genes, ]
#' ary.eff.nc <- ary.eff[!rownames(ary.eff) %in% ctrl.genes, ]
#' group.id <- substr(colnames(smp.eff.nc), 7, 7)
#'
#' redhalf.smp.eff.nc <- reduce.signal(smp.eff = smp.eff.nc,
#'                                     group.id = group.id,
#'                                     group.id.level = c("E", "V"),
#'                                     reduce.multiplier = 1/2)

"reduce.signal" <- function(smp.eff,
                                 group.id,
                                 group.id.level = c("E", "V"),
                                 reduce.multiplier = 1/2,
                                 pbset.id = NULL){
  stopifnot(nrow(smp.eff) != length(unique(rownames(smp.eff)))) # probe level
  stopifnot(length(unique(table(rownames(smp.eff)))) == 1)
  if(is.null(pbset.id)) pbset.id <- unique(rownames(smp.eff))

  n.p.u <- unique(table(rownames(smp.eff)))
  smp.eff.psl <- med.sum.pbset(smp.eff, num.per.unipbset = n.p.u)
  s.e.limma <- limma.pbset(data = smp.eff.psl,
                           group.id = group.id,
                           group.id.level = group.id.level,
                           pbset.id = pbset.id)
  de.ind <- s.e.limma$P.Value < 0.01

  sample.g1 <- rowMeans(smp.eff[rownames(smp.eff) %in% pbset.id[de.ind],
                                group.id == group.id.level[1]])
  sample.g2 <- rowMeans(smp.eff[rownames(smp.eff) %in% pbset.id[de.ind],
                                group.id == group.id.level[2]])

  half.signal <- (sample.g1 - sample.g2)*reduce.multiplier

  reduced.smp.ef.de <- cbind(smp.eff[rownames(smp.eff) %in% pbset.id[de.ind],
                                  group.id == group.id.level[1]] - half.signal,
                             smp.eff[rownames(smp.eff) %in% pbset.id[de.ind],
                                  group.id == group.id.level[2]])

  # combine and colnames, rownames back to original order
  temp <- smp.eff

  temp[rownames(temp) %in% pbset.id[de.ind], ] <- reduced.smp.ef.de[, colnames(smp.eff)]
  redhalf.sample.effect.pl.p10 <- temp; rm(temp)

  return(redhalf.sample.effect.pl.p10)
}
