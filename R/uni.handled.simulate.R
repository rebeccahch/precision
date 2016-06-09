#' Classification analysis of uniformly-handled data
#'
#' Performs classification analysis on the uniformly-handled data by reassigning samples to training and test set in Qin et al. (see reference).
#'
#' @references http://clincancerres.aacrjournals.org/content/20/13/3371.long
#' @param myseed specifies seed for random assignment using set.seed().
#' @param N number of simulation runs.
#' @param smp.eff sample effect data, rows as probes, columns as samples.
#' @param norm.list a list of strings for normalization methods compared in the simulation study;
#' built-in normalization methods includes "NN", "QN", "MN", "VSN" for "No Normalization", "Quantile Normalization", "Median Normalization", "Variance Stablizing Normalization";
#' user can provide a list of normalization methods given the functions are supplied (also see norm.funcs).
#' @param class.list a list of strings for classification methods compared in the simulation study;
#' built-in classification methods are "PAM" and "LASSO" for "prediction analysis for microarrays" and "least absolute shrinkage and selection operator";
#' user can provide a list of classification methods given the correponding model-building and predicting functions are supplied (also see class.funcs and pred.funcs).
#' @param norm.funcs a list of strings for names of user-defined normalization method functions, in the order of norm.list excluding any built-in normalization methods.
#' @param class.funcs a list of strings for names of user-defined classification model-building functions, in the order of class.list excluding any built-in classification methods.
#' @param pred.funcs a list of strings for names of user-defined classification predicting functions, in the order of class.list excluding any built-in classification methods.
#' @return benchmark analysis results with list of models built and internal and external misclassification error stored, also a list of assignment stored
#' @keywords simulation
#' @export
#' @examples
#' \dontrun{
#' smp.eff <- estimate.smp.eff(r.data = r.data.pl)
#' ctrl.genes <- unique(rownames(r.data.pl))[grep("NC", unique(rownames(r.data.pl)))]
#' smp.eff.nc <- smp.eff[!rownames(smp.eff) %in% ctrl.genes, ]
#' uni.handled.results <- uni.handled.simulate(myseed = 1, N = 3,
#'                                             smp.eff = smp.eff.nc,
#'                                             norm.list = c("NN", "QN"),
#'                                             class.list = c("PAM", "LASSO"))
#' }

"uni.handled.simulate" <- function(myseed, N, smp.eff,
                                   norm.list = c("NN", "QN"),
                                   class.list = c("PAM", "LASSO"),
                                   norm.funcs = NULL,
                                   class.funcs = NULL,
                                   pred.funcs = NULL){

  n.norm <- length(norm.list)
  n.class <- length(class.list)

  assign_store <- create.storage(design.list = "",
                                      norm.list = "",
                                      class.list = c("train", "test"),
                                      validating.sets = "assign")

  model_store <- create.storage(design.list = "",
                                     norm.list = norm.list,
                                     class.list = class.list,
                                     validating.sets = "model")

  error_store <- create.storage(design.list = "",
                                     norm.list = norm.list,
                                     class.list = class.list,
                                     validating.sets = c("internal", "external"))


  for(k in 1:N){ # each of the N simulation
    cat(k, "round seed used:", myseed + k, "\n")
    cat("- setup simulated data \n")

    #** split training and test **#
    train.ind <- sort(c(sample(which(substr(colnames(smp.eff), 7, 7) == "E"), 64),
                        sample(which(substr(colnames(smp.eff), 7, 7) == "V"), 64)))
    test.ind <- which(!colnames(smp.eff) %in% colnames(smp.eff)[train.ind])

    train <- smp.eff[, train.ind]
    test <- smp.eff[, test.ind]

    group.id.tr <- substr(colnames(train), 7, 7)
    group.id.te <- substr(colnames(test), 7, 7)

    eval(parse(text = "assign_store[[1]]$train.[[k]] <- list(train.ind)"))
    eval(parse(text = "assign_store[[1]]$test.[[k]] <- list(test.ind)"))

    cat("- preprocess data \n")
    new.norm.funcs <- switch.norm.funcs(norm.list = norm.list, norm.funcs = norm.funcs)
    for(norm.met in norm.list){
      norm.met2 <- tolower(norm.met)
      norm.func <- new.norm.funcs[norm.list == norm.met]

      if(norm.met != "NN"){
        # normalize
        eval(parse(text = paste0("temp <- ", norm.func,"(train, test)")))
        eval(parse(text = paste0("train.", norm.met2, " <- temp$train.", norm.met2)))
        eval(parse(text = paste0("test.f", norm.met2, " <- temp$test.f", norm.met2)))

        # summarize
        eval(parse(text = paste0("train.", norm.met2, ".fin <- med.sum.pbset(train.", norm.met2, ")")))
        eval(parse(text = paste0("test.f", norm.met2, ".fin <- med.sum.pbset(test.f", norm.met2, ")")))

      } else{
        eval(parse(text = "train.nn.fin <- med.sum.pbset(train)"))
        eval(parse(text = "test.fnn.fin <- med.sum.pbset(test)"))
      }

      new.class.funcs <- switch.classifier.funcs(class.list = class.list,
                                                 class.funcs = class.funcs,
                                                 pred.funcs = pred.funcs)
      for(cc in class.list){
        cat(paste0("- build ", cc, " ", norm.met, "\n"))
        class.func <- new.class.funcs$build.funcs[class.list == cc]
        pred.func <- new.class.funcs$pred.funcs[class.list == cc]

        # build
        eval(parse(text = paste0(tolower(cc), ".int.", norm.met2,
                                 " <- ", class.func,
                                 "(kfold = 5, X = train.", norm.met2,
                                 ".fin, y = group.id.tr, seed = ",
                                 myseed + k, ")")))
        # store model
        eval(parse(text = paste0("model_store['model', ][[1]][['", cc, ".",
                                 toupper(norm.met2), "']][[k]] <- list(",
                                 tolower(cc), ".int.", norm.met2, ")")))
        # store feature
        # can extract later

        # predict
        eval(parse(text = paste0(tolower(cc), ".pred.", norm.met2,
                                 " <- ", pred.func, "(", tolower(cc),
                                 ".int.", norm.met2, ", test.f", norm.met2,
                                 ".fin, group.id.te)")))
        # store coefficients
        # can extract later
        # store errors
        eval(parse(text = paste0("error_store['internal', ][[1]][['", cc, ".",
                                 toupper(norm.met2), "']][k] <- ", tolower(cc), ".int.",
                                 norm.met2, "$mc")))
        eval(parse(text = paste0("error_store['external', ][[1]][['", cc, ".",
                                 toupper(norm.met2), "']][k] <- ", tolower(cc), ".pred.",
                                 norm.met2, "$mc")))

      } # class.list loop

    } # norm.list loop

  } # simulation loop

  return(list(assign_store = assign_store,
              model_store = model_store,
              error_store = error_store))
}
