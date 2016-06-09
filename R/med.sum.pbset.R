#' Probe-set median summarization
#'
#' Summarizes probe-set using median of each unique probe, only taking in data matrix with a fixed number of probes per unique probe-set.
#'
#' @param data expression data to be summarized, rows as probes, columns as samples; rownames as probe names;
#' only accept data matrix with a fixed number of probes per unique probe-set. If the data is already on the probe-set level, no manipulation will be done.
#' @param pbset.id unique probe-set name, if not specified then use the unique probe name of the data.
#' @param num.per.unipbset number of probes for each unique probe-set; default is 10.
#' @return probe-set median summarized data
#' @importFrom stats median
#' @export
#' @keywords preprocess
#' @examples
#' r.data.psl <- med.sum.pbset(data = r.data.pl,
#'                             num.per.unipbset = 10)

"med.sum.pbset" <- function(data, pbset.id = NULL,
                            num.per.unipbset = 10) {
  stopifnot(length(unique(table(rownames(data)))) == 1)

  if(length(unique(rownames(data))) == length(rownames(data))){
    cat("Already probe-set level\n")
    data.ps <- data
  } else {
    if(is.null(pbset.id)) pbset.id <- unique(rownames(data))

    data <- data[rownames(data) %in% pbset.id, ]
    data.ps <- apply(data, 2, function(x) tapply(x, rep(1:length(unique(rownames(data))), each = num.per.unipbset), median))
    rownames(data.ps) <- pbset.id
  }
  return(data.ps)
}
