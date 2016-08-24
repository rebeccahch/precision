#' Confounding Design
#'
#' Assign arrays to samples with confounding experimental design, intentionally assigning arrays to sample groups in the order of array collection.
#' Since we split the nonuniformly-handled training dataset in a particular way so that its earlier arrays were processed by one technician and the later arrays by the other,
#' assigning the earlier arrays to one sample group and the later arrays to another
#' results in confounding handling effects with the sample groups.
#'
#' @param seed an integer used to initialize a pseudorandom number generator.
#' @param num.array number of arrays.
#' @param degree level of confounding. It must be either "complete" or "partial"
#' for complete confounding design or partial confounding design, correspondingly.
#' By default, \code{degree = "complete"}.
#' @param rev.order whether the array-to-sample-group assignment should be flipped.
#' Originally, the first half arrays are designated to be assigned to sample group 1 (the endometrial sample group)
#' and the second half to sample group 2 (the ovarian sample group).
#' If the array-to-sample-group assignment is flipped (\code{rev.order = TRUE}),
#' the first half of the array IDs will be swapped with the second half of the array IDs.
#' By default, \code{rev.order = FALSE}.
#' @return a vector of array IDs in the order of assigning to samples that are assumed to be sorted by sample group of interest
#' (first half of the samples belong to sample group 1 and second half to sample group 2).
#' As a result, the first half of the array IDs are assigned to sample group 1 and the second half of the array IDs are assigned to sample group 2.
#' @keywords study.design
#' @export
#' @examples
#'
#' # Completely confounding with reversed assignment
#' cc.rev.ind <- confounding.design(seed = 1, num.array = 128,
#'                              degree = "complete", rev.order = FALSE)
#'
#' # Partially confounding
#' pc.ind <- confounding.design(seed = 1, num.array = 128,
#'                              degree = "partial")

"confounding.design" <- function(seed, num.array,
                                 degree = "complete",
                                 rev.order = FALSE){
  stopifnot(is.numeric(seed))
  stopifnot(is.numeric(num.array))
  stopifnot(degree %in% c("complete", "partial"))
  stopifnot(num.array %% 2 == 0)

  set.seed(seed)
  if(degree == "complete"){
    g1 <- sample(1:(num.array/2))
    g2 <- sample((num.array/2 + 1):num.array)
    if(!rev.order){
      ind <- c(g1, g2)
    } else{
      ind <- c(g2, g1)
    }
  } else{ # partial
    swapsize <- ceiling(num.array/2/10)
    temp1.ind <- sample(1:(num.array/2), size = swapsize)
    temp2.ind <- sample((num.array/2+1):num.array, size = swapsize)
    g1 <- sample(1:(num.array/2))
    g2 <- sample(((num.array/2) + 1):num.array)
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
