#' Blocking Design
#'
#' Assign arrays to samples with blocking by (8-plex Agilent) array slide.
#'
#' @param seed an integer used to initialize a pseudorandom number generator.
#' @param num.smp number of arrays. It must be a multiple of 8.
#' @return a vector of array IDs in the order of assigning to samples that are assumed to be sorted by sample group of interest
#  (first half of the samples belong to group 1 and second half to group 2).
#' As a result, the first half of the array IDs are assigned to group 1 and the second half of the array IDs are assigned to group 2.
#' @export
#' @keywords study.design
#' @examples
#' blocking.design(seed = 1, num.smp = 128)
#'

"blocking.design" <- function(seed, num.smp){
  stopifnot(is.numeric(seed))
  stopifnot(is.numeric(num.smp))
  stopifnot(num.smp %% 8 == 0)

  set.seed(seed)
  g1.sample <- g2.sample <- NULL
  for(j in 1:(num.smp/8)){
    sample.number <- ((j-1)*8+1):(j*8)
    g1 <- sample(sample.number, size = 4)
    g2 <- sample.number[!sample.number %in% g1]
    g1.sample <- c(g1.sample, g1)
    g2.sample <- c(g2.sample, g2)
  }
  g1.sample <- sample(g1.sample)
  g2.sample <- sample(g2.sample)

  ind <- c(g1.sample, g2.sample)
  return(ind)
}
