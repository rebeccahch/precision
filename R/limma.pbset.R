#' Differential expression analysis of probe-set data
#'
#' Perform two-group differential expression analysis using "limma".
#'
#' @references Ritchie M., Phipson B., Wu D., Hu Y., Law C., Shi W. and Smyth G. (2015).
#' "limma powers differential expression analyses for RNA-sequencing and microarray studies."
#' Nucleic Acids Research, 43(7), pp. e47.
#' @param data dataset to be analyzed.
#' The dataset must have rows as unique probe-sets and columns as samples.
#' @param group.id a vector of sample-group labels for each sample of the dataset.
#' It must be a 2-level non-numeric factor vector.
#' @param group.id.level a vector of sample-group label level.
#' It must have two and only two elements and the first element is the reference.
#' By default, \code{group.id.level = c("E", "V")}.
#' That is in our study, we compare endometrial tumor samples to
#' ovarian tumor samples, with endometrial as our reference.
#' @param pbset.id a vector of unique probe-set names.
#' By default, \code{pbset.id = NULL} for it to be the row names of the dataset.
#' @return a data frame with differential expression analysis results,
#' group means and group standard deviations, for each unique probe-set.
#' @keywords DEA
#' @import limma
#' @importFrom stats model.matrix sd
#' @export
#' @examples
#' uhdata.psl <- med.sum.pbset(data = uhdata.pl,
#'                             num.per.unipbset = 10)
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#'
#' uhdata.psl.nc <- uhdata.psl[!rownames(uhdata.psl) %in% ctrl.genes, ]
#'
#' group.id <- substr(colnames(uhdata.psl.nc), 7, 7)
#' group.id.level <- levels(as.factor(group.id))
#'
#' limma.fit.uhdata<- limma.pbset(data = uhdata.psl.nc,
#'                                group.id = group.id,
#'                                group.id.level = group.id.level)
#'                                table(limma.fit.uhdata$P.Value < 0.01,
#'                                dnn = "DE genes")
#'

"limma.pbset" <- function(data, group.id,
                          group.id.level = c("E", "V"),
                          pbset.id = NULL){

  stopifnot(length(unique(rownames(data))) == nrow(data))
  stopifnot(is.character(group.id))
  stopifnot(group.id %in% group.id.level)

  if(is.null(pbset.id)) pbset.id <- rownames(data)

  limma.level <- factor(group.id, levels = group.id.level)
  design.mat <- model.matrix(~0 + limma.level)
  colnames(design.mat) <- group.id.level
  cont.mat <- limma::makeContrasts(contrasts = paste0(group.id.level[2], "-", group.id.level[1]),
                                   levels = design.mat)
  fit.temp <- limma::lmFit(data, design.mat)
  contr.temp <- limma::contrasts.fit(fit.temp, cont.mat)
  eb.temp <- limma::eBayes(contr.temp)
  final.temp.1 <- limma::topTable(eb.temp,number = nrow(data))
  final.temp <- final.temp.1[match(rownames(data), rownames(final.temp.1)),]

  g1.mean <- apply(data[, group.id == group.id.level[1]], 1, mean)
  g2.mean <- apply(data[, group.id == group.id.level[2]], 1, mean)
  g1.sd <- apply(data[, group.id == group.id.level[1]], 1, sd)
  g2.sd <- apply(data[, group.id == group.id.level[2]], 1, sd)

  format.t <- data.frame(final.temp, g1.mean, g2.mean, g1.sd, g2.sd)
  name.a <- ncol(final.temp) + 1
  names(format.t)[name.a:ncol(format.t)] <- c("g1.mean","g2.mean", "g1.sd", "g2.sd")

  return(format.t)
}
