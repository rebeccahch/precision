#' Stratification Design
#'
#' Assigns arrays to samples with stratification design.
#'
#' @param seed specifies seed for random assignment using set.seed().
#' @param num.smp number of samples.
#' @param batch.id sample group ID for the estimated sample effect data.
#' @return array-to-sample assignment, first half for group 1 (endometrial), second half for group 2 (ovarian)
#' @export
#' @keywords study.design
#' @examples
#' batch.id <- list(1:40, 41:64, (129:160) - 64, (161:192) - 64)
#' str.ind <- stratification.design(seed = 1, num.smp = 128,
#'                                  batch.id = batch.id)

"stratification.design" <- function(seed, num.smp,
                                    batch.id){

  stopifnot(length(unlist(batch.id)) == num.smp)
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
