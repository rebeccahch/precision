"ary.to.smp.assign.func" <- function(design, seed, n, batch.id = NULL){
  if(design == "CC+") ind <- confounding.design(seed = seed,
                                                num.smp = n,
                                                degree = "complete",
                                                rev.order = FALSE)
  if(design == "CC-") ind <- confounding.design(seed = seed,
                                                num.smp = n,
                                                degree = "complete",
                                                rev.order = TRUE)
  if(design == "PC+") ind <- confounding.design(seed = seed,
                                                num.smp = n,
                                                degree = "partial",
                                                rev.order = FALSE)
  if(design == "PC-") ind <- confounding.design(seed = seed,
                                                num.smp = n,
                                                degree = "partial",
                                                rev.order = TRUE)
  if(design == "BLK") ind <- blocking.design(seed = seed,
                                             num.smp = n)
  if(design == "STR") ind <- stratification.design(seed = seed,
                                                   num.smp = n,
                                                   batch.id = batch.id)

  return(ind)
}
