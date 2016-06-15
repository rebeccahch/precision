#' Gene type classification
#'
#' Classify genes into technical, biological or other based on the differential expression analysis results of the estimated sample and array effect data.
#'
#' @details Based on differential expression analysis, biological genes are defined as genes with \code{p-value < 0.01} on
#' all three of the overall, the training set and the test set of the sample effect data.
#' Technical genes are defined as genes with \code{p-value < 0.01} on
#' the training set of the array effect data but with \code{p-value > 0.01} on all of the following:
#' the test set of the array effect data, and the overall, the training set and the test set of the sample effect data.
#' @param smp.eff the estimated sample effect dataset. The dataset must have rows as probes and columns as samples.
#' It can only be either probe-level with a fixed number of probes per unique probe or probe-set-level.
#' @param ary.eff the estimated array effect dataset. The dataset must have rows as probes and columns as samples.
#' It must have the same dimensions and the same probe names as the estimated sample effect dataset.
#' It can only be either probe-level with a fixed number of probes per unique probe or probe-set-level.
#' @param smp.eff.train.ind a vector of sample IDs for the estimated sample effects that are assigned to the training set.
#' @param ary.eff.train.ind a vector of array IDs for the estimated array effects that are assigned to the training set.
#' @param ary.to.smp.assign a vector of indices that assign arrays to samples (see details in \code{blocking.design},
#' \code{confounding.design} or \code{stratification.design}). It must have an equal length to the number of samples
#' in the estimated sample effect dataset.
#' The first half arrays in the vector have to be assigned to the sample group 1 and the second half to sample group 2.
#' @param group.id a vector of sample-group labels for each sample of the estimated sample effect data.
#' @return a vector of gene type for each gene: -1 for technical, 0 for other, 1 for biological genes
#' @keywords classification
#' @export
#' @examples
#' \dontrun{
#' smp.eff <- estimate.smp.eff(r.data = r.data.pl)
#' ary.eff <- estimate.ary.eff(r.data = r.data.pl,
#'                             non.r.data = non.r.data.pl)
#'
#' ctrl.genes <- unique(rownames(r.data.pl))[grep("NC", unique(rownames(r.data.pl)))]
#'
#' smp.eff.nc <- smp.eff[!rownames(smp.eff) %in% ctrl.genes, ]
#' ary.eff.nc <- ary.eff[!rownames(ary.eff) %in% ctrl.genes, ]
#'
#' group.id <- substr(colnames(smp.eff.nc), 7, 7)
#'
#' smp.eff.train.ind <- colnames(smp.eff.nc)[c(sample(which(group.id == "E"), size = 64),
#' sample(which(group.id == "V"), size = 64))]
#' smp.eff.test.ind <- colnames(smp.eff.nc)[!colnames(smp.eff.nc) %in% smp.eff.train.ind]
#'
#' ary.eff.train.ind <- colnames(ary.eff.nc)[c(1:64, 129:192)]
#'
#' group.id.list <- list("all" = group.id,
#'                       "tr" = substr(smp.eff.train.ind, 7, 7),
#'                       "te" = substr(smp.eff.test.ind, 7, 7))
#'
#' ary.to.smp.assign <- list("all" = c(rep(c("E", "V"), each = 64),
#'                           rep(c("V", "E"), each = 32)),
#'                           "tr" = rep(c("E", "V"), each = 64),
#'                           "te" = rep(c("V", "E"), each = 32))
#'
#' gene.cat <- classify.gene.type(smp.eff = smp.eff.nc,
#'                                ary.eff = ary.eff.nc,
#'                                smp.eff.train.ind = smp.eff.train.ind,
#'                                ary.eff.train.ind = ary.eff.train.ind,
#'                                group.id = group.id.list,
#'                                ary.to.smp.assign = ary.to.smp.assign)
#' }

"classify.gene.type" <- function(smp.eff,
                                      ary.eff,
                                      smp.eff.train.ind,
                                      ary.eff.train.ind,
                                      group.id,
                                      ary.to.smp.assign){

  stopifnot(dim(smp.eff) == dim(ary.eff))
  stopifnot(rownames(smp.eff) == rownames(ary.eff))
  stopifnot(length(unique(table(rownames(data)))) == 1)

  if(nrow(smp.eff) != length(unique(rownames(smp.eff)))){
    cat("summarizing sample effect...\n")
    n.p.u <- unique(table(rownames(smp.eff)))
    smp.eff <- med.sum.pbset(smp.eff, num.per.unipbset = n.p.u)
  }
  if(nrow(ary.eff) != length(unique(rownames(ary.eff)))){
    cat("summarizing array effect...\n")
    n.p.u <- unique(table(rownames(ary.eff)))
    ary.eff <- med.sum.pbset(ary.eff, num.per.unipbset = n.p.u)
  }

  stopifnot(dim(smp.eff) == dim(ary.eff))

  s.e <- smp.eff
  s.e.tr <- smp.eff[, smp.eff.train.ind]
  s.e.te <- smp.eff[, !colnames(smp.eff) %in% smp.eff.train.ind]

  a.e <- ary.eff
  a.e.tr <- ary.eff[, ary.eff.train.ind]
  a.e.te <- ary.eff[, !colnames(ary.eff) %in% ary.eff.train.ind]

  s.e.limma <- limma.pbset(data = s.e, group.id = group.id$all, pbset.id = rownames(s.e))
  s.e.tr.limma <- limma.pbset(data = s.e.tr, group.id = group.id$tr, pbset.id = rownames(s.e.tr))
  s.e.te.limma <- limma.pbset(data = s.e.te, group.id = group.id$te, pbset.id = rownames(s.e.te))

  # biological genes: DE for whole sample effect dataset, training sample effect dataset and test sample effect dataset
  bio.genes <- s.e.te.limma$P.Value < 0.01 & s.e.tr.limma$P.Value < 0.01

  a.e.limma <- limma.pbset(data = a.e[, c(ary.eff.train.ind,
                                       colnames(a.e)[!colnames(a.e) %in% ary.eff.train.ind])],
                           group.id = ary.to.smp.assign$all, pbset.id = rownames(a.e))
  a.e.tr.limma <- limma.pbset(data = a.e.tr, group.id = ary.to.smp.assign$tr, pbset.id = rownames(a.e.tr))
  a.e.te.limma <- limma.pbset(data = a.e.te, group.id = ary.to.smp.assign$te, pbset.id = rownames(a.e.te))

  # technical genes: DE for training array effect dataset,
  #     non-DE for test array effect dataset
  #     non-DE for whole sample effect dataset, training sample effect dataset and test sample effect dataset
  tech.genes <- s.e.limma$P.Value > 0.01 &
    s.e.tr.limma$P.Value > 0.01 &
    s.e.te.limma$P.Value > 0.01 &
    a.e.tr.limma$P.Value < 0.01 &
    a.e.te.limma$P.Value > 0.01

  gene.cat <- ifelse(bio.genes, 1, ifelse(tech.genes, -1, 0))

  return(gene.cat)
}
