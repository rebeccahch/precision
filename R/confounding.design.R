#' Confounding Design
#'
#' Assigns arrays to samples with confounding design.
#'
#' @param seed specifies seed for random assignment using set.seed().
#' @param num.smp number of samples.
#' @param degree level of confounding; has to be either "complete" or "partial" for either "complete confounding" or "partial confounding" design; default is "complete".
#' @param rev.order FALSE indicating no reverse order, first half arrays to group 1 (endometrial), second half arrays to group 2 (ovarian); default is FALSE.
#' @return array-to-sample assignment, first half for group 1 (endometrial), second half for group 2 (ovarian)
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
