#' Blocking Design
#'
#' Assign arrays to samples with blocking experimental design by array slide.
#' For example Agilent comes with an eight-plex design.
#' That is the number of arrays in an array slide is eight;
#' within each array slide, equal number of arrays are assigned to either sample group.
#'
#' @param seed an integer used to initialize a pseudorandom number generator.
#' @param num.per.slide number of arrays per array slide. It must be a multiple of 2. By default, \code{num.per.slide = 8} for Agilent microarrays.
#' @param num.array number of arrays. It must be a multiple of \code{num.per.slide}.
#' @return a vector of array IDs in the order of assigning to samples that are assumed to be sorted by sample group of interest
#  (first half of the samples belong to sample group 1 and second half to sample group 2).
#' As a result, the first half of the array IDs are assigned to sample group 1 and the second half of the array IDs are assigned to sample group 2.
#' @export
#' @keywords study.design
#' @examples
#' blocking.design(seed = 1, num.per.slide = 8, num.array = 128)
#'

"blocking.design" <- function(seed, num.per.slide = 8, num.array){
  stopifnot(is.numeric(seed))
  stopifnot(is.numeric(num.array))
  stopifnot(num.array %% num.per.slide == 0)

  set.seed(seed)
  g1.sample <- g2.sample <- NULL
  for(j in 1:(num.array/num.per.slide)){
    sample.number <- ((j - 1)*num.per.slide + 1):(j*num.per.slide)
    g1 <- sample(sample.number, size = num.per.slide/2)
    g2 <- sample.number[!sample.number %in% g1]
    g1.sample <- c(g1.sample, g1)
    g2.sample <- c(g2.sample, g2)
  }
  g1.sample <- sample(g1.sample)
  g2.sample <- sample(g2.sample)

  ind <- c(g1.sample, g2.sample)
  return(ind)
}
