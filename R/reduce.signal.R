#' Biological signal reduction
#'
#' Reduces biological effect between sample group by a multiplier.
#'
#' @param smp.eff estimated sample effect data, rows as probes, columns as samples; can only take in probe-level data with 10 probe per unique probe.
#' @param group.id sample group ID for the estimated sample effect data.
#' @param group.id.level sample group label level; default = c("E", "V") in our studies when comparing endometrial to ovarian samples.
#' @param reduce.multiplier a multiplier specified to reduce between-sample-group signal by; default is 1/2.
#' @param pbset.id unique probe-set name, if not specified then use the unique probe name of the data.
#' @return estimated sample effect data with reduced biological signal
#' @keywords data.setup
#' @export
#' @examples
#' smp.eff <- estimate.smp.eff(r.data = r.data.pl)
#' ary.eff <- estimate.ary.eff(r.data = r.data.pl,
#'                              non.r.data = non.r.data.pl)
#' ctrl.genes <- unique(rownames(r.data.pl))[grep("NC", unique(rownames(r.data.pl)))]
#' smp.eff.nc <- smp.eff[!rownames(smp.eff) %in% ctrl.genes, ]
#' ary.eff.nc <- ary.eff[!rownames(ary.eff) %in% ctrl.genes, ]
#' group.id <- substr(colnames(smp.eff.nc), 7, 7)
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
  if(is.null(pbset.id)) pbset.id <- unique(rownames(smp.eff))

  smp.eff.psl <- med.sum.pbset(smp.eff)
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

  # sum(!rownames(temp[rownames(smp.eff) %in% pbset.id[de.ind], ]) == rownames(reduced.smp.ef.de)) # check if rownames match
  # sum(!colnames(smp.eff[rownames(smp.eff) %in% pbset.id[de.ind], ]) == colnames(reduced.smp.ef.de)) # colnames need to be matched
  # sum(!colnames(smp.eff[rownames(smp.eff) %in% pbset.id[de.ind], ]) == colnames(reduced.smp.ef.de[, colnames(smp.eff)])) # reorder like this

  temp[rownames(temp) %in% pbset.id[de.ind], ] <- reduced.smp.ef.de[, colnames(smp.eff)]
  redhalf.sample.effect.pl.p10 <- temp; rm(temp)

  return(redhalf.sample.effect.pl.p10)
}
