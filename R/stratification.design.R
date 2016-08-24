#' Stratification Design
#'
#' Assign arrays to samples with stratification,
#' a study design assigning arrays in each batch to each sample group proportionally.
#'
#' @param seed an integer used to initialize a pseudorandom number generator.
#' @param num.array number of arrays.
#' @param batch.id a list of array indices grouped by batches when data were profiled.
#' The length of the list must be equal to the number of batches in the data;
#' the number of array indices must be the same as the number of samples.
#' @return a vector of array IDs in the order of assigning to samples that are
#' assumed to be sorted by sample group of interest
#  (first half of the samples belong to group 1 and second half to group 2).
#' As a result, the first half of the array IDs are assigned to group 1
#' and the second half of the array IDs are assigned to group 2.
#' @export
#' @keywords study.design
#' @examples
#' batch.id <- list(1:40, 41:64, (129:160) - 64, (161:192) - 64)
#' str.ind <- stratification.design(seed = 1, num.array = 128,
#'                                  batch.id = batch.id)

"stratification.design" <- function(seed, num.array,
                                    batch.id){

  stopifnot(length(unlist(batch.id)) == num.array)
  set.seed(seed)

  g1.sample <- g2.sample <- NULL
  for(j in 1:length(batch.id)){
    sample.number <- batch.id[[j]]
    g1 <- sample(sample.number, size = length(sample.number)/2)
    g2 <- sample.number[!sample.number %in% g1]
    g1.sample <- c(g1.sample, g1)
    g2.sample <- c(g2.sample, g2)
  }
  g1.sample <- sample(g1.sample)
  g2.sample <- sample(g2.sample)

  ind <- c(g1.sample, g2.sample)
  return(ind)
}
