#' Classification analysis of uniformly-handled data
#'
#' Perform classification analysis on the uniformly-handled data by re-assigning samples to training and test set.
#' More details can be found in Qin et al. (see reference).
#'
#' @references Qin LX, Huang HC, Begg CB. Cautionary note on cross validation in molecular classification. Journal of Clinical Oncology. 2016.
#' @details The analysis for the uniformly-handled dataset consists of the following main steps:
#'
#' (1) randomly split the data into a training set and a test set, balanced by sample group of interest
#'
#' (2) preprocess the training data and the test data
#'
#' (3) build a classifier using the preprocessed training data
#'
#' (4) assess the mislcassification error rate of the classifier using the preprocessed test data
#'
#' This analysis is repeated for \code{N} random splits of training set and test set.
#'
#' Data preprocessing in (2) includes three steps: log2 transformation, normalization for training data
#' and frozen normalization for test data,
#' and probe-set summarization using median. Normalization methods are specified in \code{norm.list}.
#'
#' Classifier building in (3) includes choosing the tuning parameter for each method using five-fold cross-validation and
#' measuring classifier accuarcy using the misclassification error rate.
#' Classification methods are specified in \code{class.list}
#'
#' The error rate is evaluated by both external validation of test data and cross-validation of training data.
#' For user-defined normalization method or classification method, please refer to the vignette.
#'
#' @param seed an integer used to initialize a pseudorandom number generator.
#' @param N number of simulation runs.
#' @param biological.effect the estimated biological effect dataset. This dataset must have rows as probes and columns as samples.
#' @param norm.list a list of strings for normalization methods to be compared in the simulation study.
#' The built-in normalization methods includes "NN", "QN", "MN", "VSN" for "No Normalization", "Quantile Normalization",
#' "Median Normalization", "Variance Stabilizing Normalization".
#' User can provide a list of normalization methods given the functions are supplied (also see \code{norm.funcs}).
#' @param class.list a list of strings for classification methods to be compared in the simulation study.
#' The built-in classification methods are "PAM" and "LASSO" for "prediction analysis for microarrays"
#' and "least absolute shrinkage and selection operator".
#' User can provide a list of classification methods given the correponding model-building
#' and predicting functions are supplied (also see \code{class.funcs} and \code{pred.funcs}).
#' @param norm.funcs a list of strings for names of user-defined normalization method functions, in the order of \code{norm.list},
#' excluding any built-in normalization methods.
#' @param class.funcs a list of strings for names of user-defined classification model-building functions, in the order of \code{class.list},
#' excluding any built-in classification methods.
#' @param pred.funcs a list of strings for names of user-defined classification predicting functions, in the order of \code{class.list},
#' excluding any built-in classification methods.
#' @return benchmark analysis results -- a list of training-and-test-set splits, fitted models,
#' and misclassification error rates across simulation runs:
#' \item{assign_store}{random training-and-test-set splits}
#' \item{model_store}{models for each combination of normalization methods and classification methods}
#' \item{error_store}{internal and external misclassification error rates for each combination of normalization methods and classification methods}
#' @keywords simulation
#' @export
#' @examples
#' \dontrun{
#' biological.effect <- estimate.biological.effect(uhdata = uhdata.pl)
#'
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#'
#' biological.effect.nc <- biological.effect[!rownames(biological.effect) %in% ctrl.genes, ]
#'
#' uni.handled.results <- uni.handled.simulate(seed = 1, N = 3,
#'                                             biological.effect = biological.effect.nc,
#'                                             norm.list = c("NN", "QN"),
#'                                             class.list = c("PAM", "LASSO"))
#' }

"uni.handled.simulate" <- function(seed, N, biological.effect,
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
    cat(k, "round seed used:", seed + k, "\n")
    cat("- setup simulated data \n")

    #** split training and test **#
    train.ind <- sort(c(sample(which(substr(colnames(biological.effect), 7, 7) == "E"), 64),
                        sample(which(substr(colnames(biological.effect), 7, 7) == "V"), 64)))
    test.ind <- which(!colnames(biological.effect) %in% colnames(biological.effect)[train.ind])

    train <- biological.effect[, train.ind]
    test <- biological.effect[, test.ind]

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
                                 seed + k, ")")))
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
