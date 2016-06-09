#' Blocking Design
#'
#' Assigns arrays to samples with blocking design.
#'
#' @param seed specifies seed for random assignment using set.seed().
#' @param num.smp number of samples.
#' @return array-to-sample assignment, first half for group 1 (endometrial), second half for group 2 (ovarian)
#' @export
#' @keywords study.design
#' @examples
#' blocking.design(seed = 1, num.smp = 128)
#'

"blocking.design" <- function(seed, num.smp){

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
