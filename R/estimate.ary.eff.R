#' Estimated array effect
#'
#' Estimates array effect from taking the differences between the expressions of the non-randomized and the randomized data, matched by samples
#'
#' @param r.data randomized expression dataset, rows as probes, columns as samples.
#' @param non.r.data non-randomized expression dataset, rows as probes, columns as samples; must have same dimensions and same probe name as randomized data.
#' @return an estimation of the array effect
#' @keywords data.setup
#' @importFrom stats median
#' @export
#' @examples
#' ary.eff <- estimate.ary.eff(r.data = r.data.pl, non.r.data = non.r.data.pl)


"estimate.ary.eff" <- function(r.data, non.r.data){

  stopifnot(rownames(r.data) == rownames(non.r.data))
  stopifnot(dim(r.data) == dim(non.r.data))

  ## the list of 7 box 3 arrays and 2 bad arrays that require removal + median imputation
  ## from the rest of the arrays in the same slide
  rm.list <- c("JB4160V.b3","JB4387V.b3","JB4388V.b3",
               "JB4650V.b3","JB4757V.b3","JB4833V.b3",
               "GL3793V.b3","JB4933E","JB4952V")

  temp.ary.eff <- NULL
  for(i in 1:ncol(non.r.data)){
    temp.name <- substr(colnames(non.r.data)[i], 1, 7)
    b.data <- r.data[, colnames(r.data) == temp.name]
    temp.diff <- non.r.data[, i] - b.data
    temp.ary.eff <- cbind(temp.ary.eff, temp.diff)
  }
  colnames(temp.ary.eff) <- colnames(non.r.data)[1:ncol(non.r.data)]

  ary.eff <- NULL
  for(i in 1:24){
    begin.n <- (i-1)*8 + 1
    end.n <- i*8
    temp.data <- temp.ary.eff[, begin.n:end.n]
    temp.name <- colnames(temp.data)
    indi.vec <- temp.name %in% rm.list
    md.array <- apply(temp.data[, !(temp.name %in% rm.list)], 1, median)
    for(j in 1:8){
      if(indi.vec[j]){
        ary.eff <- cbind(ary.eff, md.array)
      } else{
        ary.eff <- cbind(ary.eff, temp.data[, j])
      }
    }
  }

  colnames(ary.eff) <- paste(substr(colnames(temp.ary.eff), 1, 7),
                             seq(1, ncol(non.r.data)), sep = ".")

  return(ary.eff)
}
