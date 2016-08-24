#' Probe-level data truncation to a fixed number of probes per unique probe-set
#'
#' Truncate probe-level dataset so that it has a fixed number of probes per unique probe-set.
#' We are safe to do so if the variation among replicates for the same probe is small.
#'
#' @param data probe-level dataset. The dataset must have rows as probes and columns as samples.
#' @param pbset.id a vector of unique probe-set names.
#' By default, \code{pbset.id = NULL} for it to be the row names of the dataset.
#' @param num.per.unipbset number of probes for each unique probe-set to be truncated to.
#' By default, \code{num.per.unipbset = 10}.
#' @return truncated probe-level data
#' @keywords data.setup
#' @export
#' @examples
#' uhdata.pl.p5 <- per.unipbset.truncate(data = uhdata.pl,
#' num.per.unipbset = 5)

"per.unipbset.truncate" <- function(data, pbset.id = NULL,
                                         num.per.unipbset = 10){
  if(is.null(pbset.id)) pbset.id <- unique(rownames(data))

  temp <- NULL
  for(i in pbset.id) temp <- rbind(temp, data[rownames(data) %in% i, ][1:num.per.unipbset, ])

  return(temp)
}
