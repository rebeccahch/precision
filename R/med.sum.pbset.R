#' Probe-set median summarization
#'
#' Summarize probe-set using median of each unique probe and
#' only takes in data matrix with the same number of probes per unique probe-set.
#'
#' @param data dataset to be summarized.
#' The dataset must have rows as probes and columns as samples.
#' It must be a data matrix with the same number of probes per unique probe-set.
#' If it is already on the probe-set level, no manipulation will be done.
#' @param pbset.id a vector of unique probe-set names.
#' If it is not specified, then by default it is set to be the unique probe names of the data.
#' @param num.per.unipbset number of probes for each unique probe-set.
#' By default, \code{num.per.unipbset = 10}.
#' @return probe-set median summarized data
#' @importFrom stats median
#' @export
#' @keywords preprocess
#' @examples
#' \dontrun{
#' uhdata.psl <- med.sum.pbset(data = uhdata.pl,
#'                             num.per.unipbset = 10)
#' }
#'

"med.sum.pbset" <- function(data, pbset.id = NULL,
                            num.per.unipbset = 10) {
  stopifnot(length(unique(table(rownames(data)))) == 1)

  if(length(unique(rownames(data))) == length(rownames(data))){
    cat("Already probe-set level\n")
    data.ps <- data
  } else {
    if(is.null(pbset.id)) pbset.id <- unique(rownames(data))

    data <- data[rownames(data) %in% pbset.id, ]
    data.ps <- apply(data, 2, function(x) tapply(x, rep(1:length(unique(rownames(data))),
                                                        each = num.per.unipbset), median))
    rownames(data.ps) <- pbset.id
  }
  return(data.ps)
}
