#' Classification analysis of simulation study
#'
#' Perform the simulation study in Qin et al. (see reference).
#'
#' @references Qin LX, Huang HC, Begg CB. Cautionary note on cross validation in molecular classification. Journal of Clinical Oncology. 2016
#' @details The classification anlaysis of simulation study consists of the following main steps:
#'
#' First, \code{precision.simulate} requires the training and test sets for both estimated biological effects and estimated handling effects.
#' The effects can be simulated as follows (using \code{estimate.biological.effect} and \code{estimate.handling.effect}).
#' The uniformly-handled dataset are used to approximate the biological effect for each sample,
#' and the difference between the two arrays (one from the uniformly-handled dataset and
#' the other from the nonuniformly-handled dataset, subtracting the former from the latter)
#' for the same sample are used to approximate the handling effect for each array in the nonuniformly-handled dataset.
#'
#' The samples are randomly split into a training set and a test set, balanced by tumor type (in Qin et al., training-to-test ratio is 2:1).
#' The arrays were then non-randomly split to a training set and a test set (in Qin et al., training set n = 128 -- the first 64 and last 64 arrays
#' in the order of array processing; test set n = 64 -- the middle 64 arrays).
#' This setup allows different pairings of arrays and samples by various different training-and-test-set splits.
#' Furthermore, biological signal strength and confounding level of the handling effects can be modified
#' (using \code{reduce.signal} and \code{amplify.handling.effect}).
#'
#' Second, for the training set, data are simulated through "virtual re-hybridization" (using \code{rehybridize})
#' by first assigning arrays to sample groups using a confounding design or a balanced design, and
#' then summing the biological effect for a sample and the handling effect for its assigned array.
#' Rehybridization allows us to examine the use of various array-assignment schemes, specified in \code{design.list}.
#'
#' Third, the analysis for each simulated dataset follows the same steps as described
#' for the analysis of the uniformly-handled data (also see documentation on \code{uni.handled.siumate}):
#'
#' (1) data preprocessing (normalization methods are specified in \code{norm.list} and
#' batch effects can be adjusted specified with \code{icombat}, \code{isva} and \code{iruv})
#'
#' (2) classifier training (classification methods are specified in \code{class.list})
#'
#' (3) classification error estimation using both cross-validation and external validation
#'
#' The external validation is based on the test data from the uniformly-handled dataset
#' and served as the gold standard for the misclassification error estimation.
#'
#' For a given split of samples to training set versus test set,
#' \code{N} datasets will be simulated and analyzed for each array-assignment scheme.
#' For user-defined normalization method or classification method, please refer to the vignette.
#'
#' @param seed an integer used to initialize a pseudorandom number generator.
#' @param N number of simulation runs.
#' @param biological.effect.tr the training set of the estimated biological effects. This dataset must have rows as probes and columns as samples.
#' @param biological.effect.te the test set of the estimated biological effects. This dataset must have rows as probes and columns as samples.
#' It must have the same number of probes and the same probe names as the training set of the estimated biological effects.
#' @param handling.effect.tr the training set of the estimated handling effects. This dataset must have rows as probes and columns as samples.
#' It must have the same dimensions and the same probe names as the training set of the estimated biological effects.
#' @param handling.effect.te the test set of the estimated handling effects. This dataset must have rows as probes, columns as samples.
#' It must have the same dimensions and the same probe names as the training set of the estimated handling effects.
#' @param group.id.tr a vector of sample-group labels for each sample of the training set of the estimated biological effects.
#' It must be a 2-level non-numeric factor vector.
#' @param group.id.te a vector of sample-group labels for each sample of the test set of the estimated biological effects.
#' It must be a 2-level non-numeric factor vector.
#' @param design.list a list of strings for study designs to be compared in the simulation study.
#' The built-in designs are "CC+", "CC-", "PC+", "PC-", "BLK", and "STR" for "Complete Confounding 1", "Complete Confounding 2",
#' "Partial Confounding 1", "Partial Confounding 2", "Blocking", and "Stratification" in Qin et al.
#' @param norm.list a list of strings for normalization methods to be compared in the simulation study.
#' The build-in available normalization methods are "NN", "QN", "MN", "VSN" for "No Normalization", "Quantile Normalization",
#' "Median Normalization", "Variance Stabilizing Normalization".
#' User can provide a list of normalization methods given the functions are supplied (also see norm.funcs).
#' @param class.list a list of strings for classification methods to be compared in the simulation study.
#' The built-in classification methods are "PAM" and "LASSO" for "prediction analysis for microarrays" and
#' "least absolute shrinkage and selection operator".
#' User can provide a list of classification methods given the correponding model-building and
#' predicting functions are supplied (also see \code{class.funcs} and \code{pred.funcs}).
#' @param batch.id a list of array indices grouped by batches when data were profiled.
#' The length of the list must be equal to the number of batches in the data;
#' the number of array indices must be the same as the number of samples.
#' This is required if stratification study design is specified in \code{design.list}; otherwise \code{batch.id = NULL}.
#' @param icombat an indicator for combat adjustment. By default, \code{icombat = FALSE} for no ComBat adjustment.
#' @param isva an indicator for sva adjustment. By default, \code{isva = FALSE} for no sva adjustment.
#' @param iruv an indicator for RUV-4 adjustment. By default, \code{iruv = FALSE} for no RUV-4 adjustment.
#' @param biological.effect.tr.ctrl the training set of the negative-control probe biological effect data if \code{iruv = TRUE}.
#' This dataset must have rows as probes and columns as samples.
#' It also must have the same number of samples and the same sample names as \code{biological.effect.tr}.
#' @param handling.effect.tr.ctrl the training set of the negative-control probe handling effect data if \code{iruv = TRUE}.
#' This dataset must have rows as probes and columns as samples.
#' It also must have the same dimensions and the same probe names as \code{biological.effect.tr.ctrl}.
#' @param norm.funcs a list of strings for names of user-defined normalization method functions, in the order of \code{norm.list},
#' excluding any built-in normalization methods.
#' @param class.funcs a list of strings for names of user-defined classification model-building functions, in the order of \code{class.list},
#' excluding any built-in classification methods.
#' @param pred.funcs a list of strings for names of user-defined classification predicting functions, in the order of \code{class.list},
#' excluding any built-in classification methods.
#' @return simulation study results -- a list of array-to-sample assignments, fitted models,
#' and misclassification error rates across simulation runs:
#' \item{assign_store}{array-to-sample assignments for each study design}
#' \item{model_store}{models for each combination of study designs, normalization methods, and classification methods}
#' \item{error_store}{internal and external misclassification error rates for each combination of study designs,
#' normalization methods, and classification methods}
#' @keywords simulation
#' @export
#' @examples
#' \dontrun{
#' set.seed(101)
#' biological.effect <- estimate.biological.effect(uhdata = uhdata.pl)
#' handling.effect <- estimate.handling.effect(uhdata = uhdata.pl,
#'                              nuhdata = nuhdata.pl)
#'
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#'
#' biological.effect.nc <- biological.effect[!rownames(biological.effect) %in% ctrl.genes, ]
#' handling.effect.nc <- handling.effect[!rownames(handling.effect) %in% ctrl.genes, ]
#'
#' group.id <- substr(colnames(biological.effect.nc), 7, 7)
#'
#' # randomly split biological effect data into training and test set with
#' # equal number of endometrial and ovarian samples
#' biological.effect.train.ind <- colnames(biological.effect.nc)[c(sample(which(group.id == "E"), size = 64),
#'                                           sample(which(group.id == "V"), size = 64))]
#' biological.effect.test.ind <- colnames(biological.effect.nc)[!colnames(biological.effect.nc) %in% biological.effect.train.ind]
#' biological.effect.train.test.split =
#'   list("tr" = biological.effect.train.ind,
#'        "te" = biological.effect.test.ind)
#'
#' # non-randomly split handling effect data into training and test set
#' handling.effect.train.test.split =
#'   list("tr" = c(1:64, 129:192),
#'        "te" = 65:128)
#'
#' biological.effect.nc.tr <- biological.effect.nc[, biological.effect.train.ind]
#' biological.effect.nc.te <- biological.effect.nc[, biological.effect.test.ind]
#' handling.effect.nc.tr <- handling.effect.nc[, c(1:64, 129:192)]
#' handling.effect.nc.te <- handling.effect.nc[, 65:128]
#'
#' # Simulation without batch adjustment
#' precision.results <- precision.simulate(seed = 1, N = 3,
#'                                         biological.effect.tr = biological.effect.nc.tr,
#'                                         biological.effect.te = biological.effect.nc.te,
#'                                         handling.effect.tr = handling.effect.nc.tr,
#'                                         handling.effect.te = handling.effect.nc.te,
#'                                         group.id.tr = substr(colnames(biological.effect.nc.tr), 7, 7),
#'                                         group.id.te = substr(colnames(biological.effect.nc.te), 7, 7),
#'                                         design.list = c("PC-", "STR"),
#'                                         norm.list = c("NN", "QN"),
#'                                         class.list = c("PAM", "LASSO"),
#'                                         batch.id = list(1:40,
#'                                                         41:64,
#'                                                         (129:160) - 64,
#'                                                         (161:192) - 64))
#'
#' # Simulation with RUV-4 batch adjustment
#' biological.effect.ctrl <- biological.effect[rownames(biological.effect) %in% ctrl.genes, ]
#' handling.effect.ctrl <- handling.effect[rownames(handling.effect) %in% ctrl.genes, ]
#'
#' biological.effect.tr.ctrl <- biological.effect.ctrl[, biological.effect.train.test.split$tr]
#' handling.effect.tr.ctrl <- handling.effect.ctrl[, handling.effect.train.test.split$tr]
#'
#' precision.ruv4.results <- precision.simulate(seed = 1, N = 3,
#'                                              biological.effect.tr = biological.effect.nc.tr,
#'                                              biological.effect.te = biological.effect.nc.te,
#'                                              handling.effect.tr = handling.effect.nc.tr,
#'                                              handling.effect.te = handling.effect.nc.te,
#'                                              group.id.tr = substr(colnames(biological.effect.nc.tr), 7, 7),
#'                                              group.id.te = substr(colnames(biological.effect.nc.te), 7, 7),
#'                                              design.list = c("PC-", "STR"),
#'                                              norm.list = c("NN", "QN"),
#'                                              class.list = c("PAM", "LASSO"),
#'                                              batch.id = list(1:40,
#'                                                              41:64,
#'                                                              (129:160) - 64,
#'                                                              (161:192) - 64),
#'                                              iruv = TRUE,
#'                                              biological.effect.tr.ctrl = biological.effect.tr.ctrl,
#'                                              handling.effect.tr.ctrl = handling.effect.tr.ctrl)
#' }
#'
"precision.simulate" <- function(seed, N,
                                 biological.effect.tr, biological.effect.te,
                                 handling.effect.tr, handling.effect.te,
                                 group.id.tr, group.id.te,
                                 design.list = c("CC+", "CC-", "PC+", "PC-"),
                                 norm.list = c("NN", "QN"),
                                 class.list = c("PAM", "LASSO"),
                                 batch.id = NULL,
                                 icombat = FALSE, isva = FALSE, iruv = FALSE,
                                 biological.effect.tr.ctrl = NULL,
                                 handling.effect.tr.ctrl = NULL,
                                 norm.funcs = NULL,
                                 class.funcs = NULL,
                                 pred.funcs = NULL){
  stopifnot(design.list %in% c("CC+", "CC-", "PC+", "PC-", "BLK", "STR"))

  n.design <- length(design.list)
  n.norm <- length(norm.list)
  n.class <- length(class.list)

  assign_store <- create.storage(design.list = design.list,
                                 norm.list = "",
                                 class.list = "train",
                                 validating.sets = "assign")

  model_store <- create.storage(design.list = design.list,
                                norm.list = norm.list,
                                class.list = class.list,
                                validating.sets = "model")

  error_store <- create.storage(design.list = design.list,
                                norm.list = norm.list,
                                class.list = class.list,
                                validating.sets = c("internal", "external"))


  for(k in 1:N){ # each of the N simulation
    cat(k, "round seed used:", seed + k, "\n")
    #cat("- setup simulated data \n")
    for(dd in design.list){
      if(dd == "CC+") cc1.ind <- confounding.design(seed = seed + k, num.array = ncol(handling.effect.tr), degree = "complete", rev.order = FALSE)
      if(dd == "CC-") cc2.ind <- confounding.design(seed = seed + k, num.array = ncol(handling.effect.tr), degree = "complete", rev.order = TRUE)
      if(dd == "PC+") pc1.ind <- confounding.design(seed = seed + k, num.array = ncol(handling.effect.tr), degree = "partial", rev.order = FALSE)
      if(dd == "PC-") pc2.ind <- confounding.design(seed = seed + k, num.array = ncol(handling.effect.tr), degree = "partial", rev.order = TRUE)
      if(dd == "BLK") blk.ind <- blocking.design(seed = seed + k, num.array = ncol(handling.effect.tr))
      if(dd == "STR") str.ind <- stratification.design(seed = seed + k, num.array = ncol(handling.effect.tr), batch.id = batch.id)

      dd2 <- tolower(gsub("[-]", "2", gsub("[+]", "1", dd)))

      param <- ifelse(icombat, ", icombat = TRUE",
                      ifelse(isva, ", isva = TRUE",
                             ifelse(iruv, ", iruv = TRUE,
                                    biological.effect.ctrl = biological.effect.tr.ctrl,
                                    handling.effect.ctrl = handling.effect.tr.ctrl", "")))

      eval(parse(text = paste0(dd2, ".tr <- rehybridize(biological.effect = biological.effect.tr,
                               handling.effect = handling.effect.tr, group.id = group.id.tr,
                               array.to.sample.assign = ", dd2, ".ind", param, ")")))

      eval(parse(text = paste0("assign_store['assign', ][['", dd,
                               "']]$train.[[k]] <- list(", dd2, ".ind)")))

      if(isva){ # isva has special output from rehybridize()

        eval(parse(text = paste0("fsvaobj <- fsva(", dd2, ".tr$trainData, ",
                                 dd2, ".tr$trainMod, ",
                                 dd2, ".tr$trainSV, biological.effect.te)")))

        eval(parse(text = paste0(dd2, ".tr <- fsvaobj$db # train.sva")))
        eval(parse(text = "biological.effect.te2 <- fsvaobj$new # test.fsva"))
      } else {
        biological.effect.te2 <- biological.effect.te
      }

      #cat("- preprocess data \n")
      new.norm.funcs <- switch.norm.funcs(norm.list = norm.list, norm.funcs = norm.funcs)
      for(norm.met in norm.list){
        norm.met2 <- tolower(norm.met)
        norm.func <- new.norm.funcs[norm.list == norm.met]

        if(norm.met != "NN"){
          # normalize
          eval(parse(text = paste0("temp <- ", norm.func,"(", dd2, ".tr, biological.effect.te2)")))
          eval(parse(text = paste0(dd2, ".tr.", norm.met2, " <- temp$train.", norm.met2)))
          eval(parse(text = paste0("gs.", dd2, ".f", norm.met2, " <- temp$test.f", norm.met2)))

          # summarize
          eval(parse(text = paste0(dd2, ".tr.", norm.met2, ".fin <- med.sum.pbset(", dd2, ".tr.", norm.met2, ")")))
          eval(parse(text = paste0("gs.", dd2, ".f", norm.met2, ".fin <- med.sum.pbset(gs.", dd2, ".f", norm.met2, ")")))

        } else{
          eval(parse(text = paste0(dd2, ".tr.nn.fin <- med.sum.pbset(", dd2, ".tr)")))
          #eval(parse(text = paste0("gs.", dd2, ".fnn.fin <- med.sum.pbset(biological.effect.te2)")))
          eval(parse(text = paste0("temp <- quant.norm(", dd2, ".tr, biological.effect.te2)")))
          eval(parse(text = paste0("gs.", dd2, ".fqn <- temp$test.fqn")))
          eval(parse(text = paste0("gs.", dd2, ".fqn.fin <- med.sum.pbset(gs.", dd2, ".fqn)")))
        }

        new.class.funcs <- switch.classifier.funcs(class.list = class.list,
                                                   class.funcs = class.funcs,
                                                   pred.funcs = pred.funcs)

        for(cc in class.list){
          cat(paste0("- build ", cc, " ", dd, " & ", norm.met, "\n"))
          class.func <- new.class.funcs$build.funcs[class.list == cc]
          pred.func <- new.class.funcs$pred.funcs[class.list == cc]

          # build
          eval(parse(text = paste0(dd2, ".", tolower(cc), ".int.", norm.met2,
                                   " <- ", class.func,
                                   "(kfold = 5, X = ", dd2, ".tr.", norm.met2,
                                   ".fin, y = group.id.tr,
                                   seed = ", seed + k, ")")))
          # store model
          eval(parse(text = paste0("model_store['model', ][['", dd, "']][['", cc, ".",
                                   toupper(norm.met2), "']][[k]] <- list(", dd2, ".",
                                   tolower(cc), ".int.", norm.met2, ")")))
          # store feature
          # can extract later

          # predict
          eval(parse(text = paste0(dd2, ".", tolower(cc), ".pred.", norm.met2,
                                   " <- ", pred.func, "(", dd2, ".", tolower(cc),
                                   ".int.", norm.met2, ", gs.", dd2,
                                   ".fqn.fin, group.id.te)")))
          # store coefficients
          # can extract later
          # store errors
          eval(parse(text = paste0("error_store['internal', ][['", dd, "']][['", cc, ".",
                                   toupper(norm.met2), "']][k] <- ", dd2, ".",
                                   tolower(cc), ".int.", norm.met2, "$mc")))
          eval(parse(text = paste0("error_store['external', ][['", dd, "']][['", cc, ".",
                                   toupper(norm.met2), "']][k] <- ", dd2, ".",
                                   tolower(cc), ".pred.", norm.met2, "$mc")))

        } # class.list loop

      } # norm.list loop

    } # design.list loop

  } # simulation loop

  return(list(assign_store = assign_store,
              model_store = model_store,
              error_store = error_store))
}
