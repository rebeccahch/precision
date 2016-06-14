#' Confounding Design
#'
#' Assign arrays to samples with confounding design, intentionally assigning arrays to sample groups in the order of array collection.
#' Since the non-uniformly-handled data had the earlier arrays processed by one technician and the later arrays processed by another,
#' assigning the earlier arrays to one sample group and the later arrays to another
#' results in confounding handling effects with the sample groups.
#'
#' @param seed an integer used to initialize a pseudorandom number generator.
#' @param num.smp number of arrays.
#' @param degree level of confounding. It must be either "complete" or "partial"
#' for "complete confounding" or "partial confounding" design.
#' By default, \code{degree = "complete"}.
#' @param rev.order whether the array-to-sample-group assignment should be flipped.
#' Originally the first half arrays are designated to be assigned to group 1 (endometrial sample group)
#' and the second half to group 2 (ovarian sample group).
#' If the array-to-sample-group assignment is flipped (rev.order = TRUE),
#' the first half of the array IDs will be swapped with the second half of the array IDs.
#' By default, \code{rev.order = FALSE}.
#' @return a vector of array IDs in the order of assigning to samples that are assumed to be sorted by sample group of interest
#  (first half of the samples belong to group 1 and second half to group 2).
#' As a result, the first half of the array IDs are assigned to group 1 and the second half of the array IDs are assigned to group 2.
#' @keywords study.design
#' @export
#' @examples
#' cc.ind <- confounding.design(seed = 1, num.smp = 128,
#'                              degree = "complete", rev.order = FALSE)
#' cc.ind <- confounding.design(seed = 1, num.smp = 128,
#'                              degree = "complete", rev.order = FALSE)

"confounding.design" <- function(seed, num.smp,
                                 degree = "complete",
                                 rev.order = FALSE){
  stopifnot(is.numeric(seed))
  stopifnot(is.numeric(num.smp))
  stopifnot(degree %in% c("complete", "partial"))
  stopifnot(num.smp %% 2 == 0)

  set.seed(seed)
  if(degree == "complete"){
    g1 <- sample(1:(num.smp/2))
    g2 <- sample((num.smp/2 + 1):num.smp)
    if(!rev.order){
      ind <- c(g1, g2)
    } else{
      ind <- c(g2, g1)
    }
  } else{ # partial
    swapsize <- ceiling(num.smp/2/10)
    temp1.ind <- sample(1:(num.smp/2), size = swapsize)
    temp2.ind <- sample((num.smp/2+1):num.smp, size = swapsize)
    g1 <- sample(1:(num.smp/2))
    g2 <- sample(((num.smp/2) + 1):num.smp)
    g1[g1 %in% temp1.ind] <- temp2.ind
    g2[g2 %in% temp2.ind] <- temp1.ind
    rm(temp1.ind, temp2.ind)

    if(!rev.order){
      ind <- c(g1, g2)
    } else{
      ind <- c(g2, g1)
    }
  }
  return(ind)
}
