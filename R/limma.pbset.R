#' Differential expression analysis of probe-set data
#'
#' Performs two-group differential expression analysis using "limma".
#'
#' @param data expression dataset to be differentially expression analyzed, rows as unique probe-sets, columns as samples.
#' @param group.id sample group label; must be a 2-level non-numeric factor vector.
#' @param group.id.level sample group label level, the first one being the reference level; default = c("E", "V") in our studies when comparing endometrial to ovarian samples.
#' @param pbset.id unique probe-set name; default is NULL, the rownames of the dataset.
#' @return differential expression anlysis results, group means, group standard deviations
#' @keywords DEA
#' @import limma
#' @importFrom stats model.matrix sd
#' @export
#' @examples
#' r.data.psl <- med.sum.pbset(data = r.data.pl,
#'                             num.per.unipbset = 10)
#' ctrl.genes <- unique(rownames(r.data.pl))[grep("NC", unique(rownames(r.data.pl)))]
#' r.data.psl.nc <- r.data.psl[!rownames(r.data.psl) %in% ctrl.genes, ]
#' group.id <- substr(colnames(r.data.psl.nc), 7, 7)
#' group.id.level <- levels(as.factor(group.id))
#' limma.fit.r.data<- limma.pbset(data = r.data.psl.nc,
#'                                group.id = group.id,
#'                                group.id.level = group.id.level)
#'                                table(limma.fit.r.data$P.Value < 0.01, dnn = "DE genes")

"limma.pbset" <- function(data, group.id,
                          group.id.level = c("E", "V"),
                          pbset.id = NULL){

  stopifnot(length(unique(rownames(data))) == nrow(data))
  stopifnot(is.character(group.id))
  stopifnot(group.id %in% group.id.level)

  if(is.null(pbset.id)) pbset.id <- rownames(data)

  limma.level <- factor(group.id,levels = group.id.level)
  design.mat <- model.matrix(~0 + limma.level)
  colnames(design.mat) <- group.id.level
  cont.mat <- limma::makeContrasts(contrasts = paste0(group.id.level[2], "-", group.id.level[1]),
                                   levels = design.mat)
  fit.temp <- limma::lmFit(data, design.mat)
  contr.temp <- limma::contrasts.fit(fit.temp, cont.mat)
  eb.temp <- limma::eBayes(contr.temp)
  final.temp.1 <- limma::topTable(eb.temp,number = nrow(data))
  final.temp <- final.temp.1[match(rownames(data), rownames(final.temp.1)),]
  ## format and organize the limma result
  g1.mean <- apply(data[, group.id == group.id.level[1]], 1, mean)
  g2.mean <- apply(data[, group.id == group.id.level[2]], 1, mean)
  g1.sd <- apply(data[, group.id == group.id.level[1]], 1, sd)
  g2.sd <- apply(data[, group.id == group.id.level[2]], 1, sd)

  format.t <- data.frame(final.temp, g1.mean, g2.mean, g1.sd, g2.sd)# ,q.v)
  name.a <- ncol(final.temp) + 1
  names(format.t)[name.a:ncol(format.t)] <- c("g1.mean","g2.mean", "g1.sd", "g2.sd")

  return(format.t)
}
