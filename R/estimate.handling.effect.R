#' Estimated handling effects
#'
#' Estimate handling effects of an array of the nonuniformly-handled dataset
#' by taking the difference between its data and the data of
#' its matched array in the uniformly-handled dataset.
#'
#' @param uhdata the uniformly-handled dataset.
#' The dataset must have rows as probes and columns as samples.
#' @param nuhdata the nonuniformly-handled dataset.
#' The dataset must have rows as probes and columns as samples and the same dimensions and
#' the same probe names as the uniformly-handled dataset.
#' @return an estimation of the handling effects
#' @keywords data.setup
#' @importFrom stats median
#' @export
#' @examples
#' handling.effect <- estimate.handling.effect(uhdata = uhdata.pl, nuhdata = nuhdata.pl)
#'


"estimate.handling.effect" <- function(uhdata, nuhdata){

  stopifnot(rownames(uhdata) == rownames(nuhdata))
  stopifnot(dim(uhdata) == dim(nuhdata))

  ## the list of 7 box 3 arrays and 2 bad arrays that require removal + median imputation
  ## from the rest of the arrays in the same slide
  rm.list <- c("JB4160V.b3","JB4387V.b3","JB4388V.b3",
               "JB4650V.b3","JB4757V.b3","JB4833V.b3",
               "GL3793V.b3","JB4933E","JB4952V")

  temp.handling.effect <- NULL
  for(i in 1:ncol(nuhdata)){
    temp.name <- substr(colnames(nuhdata)[i], 1, 7)
    temp.data <- uhdata[, colnames(uhdata) == temp.name]
    temp.diff <- nuhdata[, i] - temp.data
    temp.handling.effect <- cbind(temp.handling.effect, temp.diff)
  }
  colnames(temp.handling.effect) <- colnames(nuhdata)[1:ncol(nuhdata)]

  handling.effect <- NULL
  for(i in 1:24){
    begin.n <- (i-1)*8 + 1
    end.n <- i*8
    temp.data <- temp.handling.effect[, begin.n:end.n]
    temp.name <- colnames(temp.data)
    indi.vec <- temp.name %in% rm.list
    md.array <- apply(temp.data[, !(temp.name %in% rm.list)], 1, median)
    for(j in 1:8){
      if(indi.vec[j]){
        handling.effect <- cbind(handling.effect, md.array)
      } else{
        handling.effect <- cbind(handling.effect, temp.data[, j])
      }
    }
  }

  colnames(handling.effect) <- paste(substr(colnames(temp.handling.effect), 1, 7),
                             seq(1, ncol(nuhdata)), sep = ".")

  return(handling.effect)
}
