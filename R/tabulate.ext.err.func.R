"tabulate.ext.err.func" <- function(pred.obj, obs.grp)
  return(1 - sum(diag(table(pred.obj, obs.grp)))/length(obs.grp))
