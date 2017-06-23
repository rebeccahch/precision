#' Plot misclassification error rates from PRECISION (both non-FLEX and FLEX) output
#'
#' Plot error rates from PRECISION output that is being extracted by \code{extract.precision.error}
#'
#' @param data a resulted PRECISION output from \code{extract.precision.error}.
#' @param design.order indices in the order of what user would like to plot from the resulted output from \code{extract.precision.error}. By default, every simulation is plotted.
#' @param mytitle the title of the plot.
#' @param ylim y limit. By default, \code{ylim = c(0, 0.5)}.
#' @param iflex whether data if a PRECISION output from \code{precision.simulate} or \code{precision.simulate.flex}. By default, \code{iflex = TRUE} indicating that the data are from \code{precision.simulate.flex}.
#' @return a boxplot of misclassification error rates
#' @keywords result.manage
#' @export
#' @examples
#' \dontrun{
#' ## PRECISION output
#' plot.precision(data = precision.results.err.df, 
#'   design.order = 1:10,
#'   mytitle = "My PRECISION results", iflex = FALSE)
#' 
#' ## PRECISION FLEX output
#' plot.precision(data = precision.results.flex.err.df, 
#'   design.order = 1:10,
#'   mytitle = "My PRECISION FLEX results")
#'           
#' # exclude ext.uh
#' plot.precision(data = precision.results.flex.err.df[!names(precision.results.flex.err.df) %in% "ext.uh"], 
#'   design.order = 1:10,
#'   mytitle = "My PRECISION FLEX results")
#' }
#'
"plot.precision" <- function(data,
                             mytitle = "PRECISION results (misclass. error rates)",
                             design.order = NULL,
                             ylim = c(0, 0.5), iflex = TRUE){
  
  if(is.null(design.order)) design.order <- 1:ncol(data[[1]])
  
  # set up frame
  boxplot(data[[1]][, design.order],
          at = 1:length(design.order) - 0.4,
          las = 2, cex.axis = 0.8,
          ylab = "Misclassification error rate",
          xlim = c(0.4, length(design.order) + 0.2),
          ylim = ylim, boxwex = 0.35,
          col = NULL, border = NA, xaxt = "n", pch = 4,
          main = mytitle)
  
  if(iflex){
    if(!is.null(data$ext.uh)){
      boxplot(data$ext.uh[, design.order], add = TRUE, 
              at = 1:length(design.order) - 0.4,
              las = 2, cex.axis = 0.8,
              boxwex = 0.35,
              col = adjustcolor("skyblue", alpha.f = 0.6),
              border = "deepskyblue3", xaxt = "n", pch = 4)
    }
    
    if(!is.null(data$ext.sim.nuh)){
      boxplot(data$ext.sim.nuh[, design.order],
              at = 1:length(design.order) - 0.4,
              las = 2, add = TRUE, border = "darkgrey", boxwex = 0.20,
              col = adjustcolor("lightgrey", alpha.f = 0.6),
              xaxt = "n", yaxt = "n", pch = 3)
    }
    
    if(!is.null(data$int)){
      boxplot(data$int[, design.order],
              at = 1:ncol(data$int[, design.order]),
              las = 2, add = TRUE, border = "red", col = NA, boxwex = 0.35,
              xaxt = "n", yaxt = "n", pch = 2)
    }
  } else{
    
    if(!is.null(data$external)){
      boxplot(data$external[, design.order], add = TRUE, 
              at = 1:length(design.order) - 0.4,
              las = 2, cex.axis = 0.8,
              boxwex = 0.35,
              col = adjustcolor("skyblue", alpha.f = 0.6),
              border = "deepskyblue3", xaxt = "n", pch = 4)
    }
    
    if(!is.null(data$internal)){
      boxplot(data$internal[, design.order],
              at = 1:ncol(data$internal[, design.order]),
              las = 2, add = TRUE, border = "red", col = NA, boxwex = 0.35,
              xaxt = "n", yaxt = "n", pch = 2)
    }
  }
  
  grid(nx = NA, ny = NULL, lwd = 2)
  labs <- names(data[[1]][, design.order])
  axis(1, cex.axis = 0.6, at = 1:length(labs) - 0.2,
       col = NA, las = 2,
       labels = labs)
  abline(v = 1:(ncol(data[[1]][, design.order]) - 1) + 0.3,
         lty = "dashed", col = "dimgrey")
  
  legend.names <- as.character(sapply(names(data), 
    function(x) switch(x, 
      "int" = "Cross Validation", 
      "internal" = "Cross Validation", 
      "ext.uh" = "E.V. with UH", 
      "external" = "E.V. with UH", 
      "ext.sim.nuh" = "E.V. with Non-UH")))
    
  border.cols <- as.character(sapply(names(data), 
    function(x) switch(x, 
      "int" = "red", 
      "internal" = "red", 
      "ext.uh" = "deepskyblue3", 
      "external" = "deepskyblue3", 
      "ext.sim.nuh" = "darkgrey")))
  
  fill.cols <- as.character(sapply(names(data), 
    function(x) switch(x, 
      "int" = NA, 
      "internal" = NA, 
      "ext.uh" = adjustcolor("skyblue", alpha.f = 0.6), 
      "external" = adjustcolor("skyblue", alpha.f = 0.6), 
      "ext.sim.nuh"  = adjustcolor("lightgrey", alpha.f = 0.6))))
    
  legend("topleft", cex = 0.5,
         legend = legend.names,
         border = border.cols,
         fill = fill.cols)  
  
}

