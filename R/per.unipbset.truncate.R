#' Classification analysis of uniformly-handled data
#'
#' Performs classification analysis on the uniformly-handled data by reassigning samples to training and test set in Qin et al. (see reference).
#'
#' @references http://clincancerres.aacrjournals.org/content/20/13/3371.long
#' @param data expression data, rows as probes, columns as samples.
#' @param pbset.id unique probe-set name; default is NULL, the rownames of the dataset.
#' @param num.per.unipbset number of probes for each unique probe-set; default is 10.
#' @return benchmark analysis results with list of models built and internal and external misclassification error stored, also a list of assignment stored
#' @keywords data.setup
#' @export
#' @examples
#' r.data.pl.p5 <- per.unipbset.truncate(data = r.data.pl,
#' num.per.unipbset = 5)

"per.unipbset.truncate" <- function(data, pbset.id = NULL,
                                         num.per.unipbset = 10){
  if(is.null(pbset.id)) pbset.id <- unique(rownames(data))

  temp <- NULL
  for(i in pbset.id) temp <- rbind(temp, data[rownames(data) %in% i, ][1:num.per.unipbset, ])

  return(temp)
}
