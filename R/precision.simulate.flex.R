#' Classification analysis of simulation study (with more flexibility)
#'
#' Perform the simulation study similar to the one in Qin et al., but allow different combinations of study designs and normalization methods on training set and test set and support one internal validation and two external validations - one using uniformly-handled test set and the other one using nonuniformly-handled test set.
#'
#' @references Qin LX, Huang HC, Begg CB. Cautionary note on cross validation in molecular classification. Journal of Clinical Oncology. 2016
#' @details The main steps of the classification anlaysis of simulation study are explained in \code{precision.simulate} in details. This function includes more flexible functionalities such as allowing different combinations of study designs and normalization methods on training and test sets. For instance, user can now compare the effect of frozen quantile normalization by running two simulations: 1) quantile normalization on both training set and test set and 2) quantile normalization on training set but forzen quantile normalization on test set. Or user can compare the effect of forzen normaliation on test set by varying frozen normalization methods on test set but fixing normalization method on traing set for two simulations. Another functionality is that user can include a third external validation method - using the nonuniformly-handled test set as the external independent set. The nonuniformly-handled test set is simulated the same way as the nonuniformly-handled training set, using user-specified study design(s).
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
#' @param design.tr.list a list of strings for study designs on the training set to be compared in the simulation study.
#' The built-in designs are "CC+", "CC-", "PC+", "PC-", "BLK", and "STR" for "Complete Confounding 1", "Complete Confounding 2",
#' "Partial Confounding 1", "Partial Confounding 2", "Blocking", and "Stratification" in Qin et al.
#' @param design.te.list a list of strings for study designs on the test set to be compared in the simulation study.
#' It must have the same length as \code{design.tr.list}. See \code{design.tr.list} for the built-in designs.
#' @param norm.tr.list a list of strings for normalization methods on the training set to be compared in the simulation study.
#' It must have the same length as \code{design.tr.list} and \code{design.te.list}.
#' The build-in available normalization methods are "NN", "QN", "MN", "VSN" for "No Normalization", "Quantile Normalization",
#' "Median Normalization", "Variance Stabilizing Normalization".
#' User can provide a list of normalization methods given the functions are supplied (also see \code{norm.tr.funcs}).
#' @param norm.te.list a list of strings for normalization methods on the test set to be compared in the simulation study.
#' It must have the same length as \code{norm.te.list}.
#' See \code{norm.tr.list} for the build-in available normalization methods.
#' User can provide a list of normalization methods given the functions are supplied (also see \code{norm.tr.funcs}).
#' @param class.list a list of strings for classification methods to be compared in the simulation study.
#' The built-in classification methods are "PAM" and "LASSO" for "prediction analysis for microarrays" and
#' "least absolute shrinkage and selection operator".
#' User can provide a list of classification methods given the correponding model-building and
#' predicting functions are supplied (also see \code{class.funcs} and \code{pred.funcs}).
#' @param valid.list a list of strings for validation methods to be compared in the simulation study. The built-in validation methods are:
#' \code{int}, \code{ext.uh}, and \code{ext.sim.nuh} which respectively represent
#' internal validation, external validation using uniformly-handled test set
#' (i.e., with biological effects only), and external validation
#' using nonuniformly-handled test set.
#' By default, \code{valid.list = c("int", "ext.uh", "ext.sim.nuh")}.
#' @param batch.id.tr a list of array indices grouped by batches when training data were profiled.
#' The length of the list must be equal to the number of batches in the training data;
#' the number of array indices must be the same as the number of samples.
#' This is required if stratification study design is specified in \code{design.tr.list}; otherwise \code{batch.id.tr = NULL}.
#' @param batch.id.te a list of array indices grouped by batches when test data were profiled.
#' The length of the list must be equal to the number of batches in the test data;
#' the number of array indices must be the same as the number of samples.
#' This is required if stratification study design is specified in \code{design.te.list}; otherwise \code{batch.id.te = NULL}.
#' @param icombat an indicator for combat adjustment. By default, \code{icombat = FALSE} for no ComBat adjustment.
#' @param isva an indicator for sva adjustment. By default, \code{isva = FALSE} for no sva adjustment.
#' @param iruv an indicator for RUV-4 adjustment. By default, \code{iruv = FALSE} for no RUV-4 adjustment.
#' @param biological.effect.tr.ctrl the training set of the negative-control probe biological effect data if \code{iruv = TRUE}.
#' This dataset must have rows as probes and columns as samples.
#' It also must have the same number of samples and the same sample names as \code{biological.effect.tr}.
#' @param handling.effect.tr.ctrl the training set of the negative-control probe handling effect data if \code{iruv = TRUE}.
#' This dataset must have rows as probes and columns as samples.
#' It also must have the same dimensions and the same probe names as \code{biological.effect.tr.ctrl}.
#' @param norm.tr.funcs a list of strings for names of user-defined normalization method functions for the training set, in the order of \code{norm.tr.list},
#' excluding any built-in normalization methods.
#' @param norm.te.funcs a list of strings for names of user-defined normalization method functions for the test set, in the order of \code{norm.te.list},
#' excluding any built-in normalization methods.
#' @param class.funcs a list of strings for names of user-defined classification model-building functions, in the order of \code{class.list},
#' excluding any built-in classification methods.
#' @param pred.funcs a list of strings for names of user-defined classification predicting functions, in the order of \code{class.list},
#' excluding any built-in classification methods.
#' @return simulation study results -- a list of array-to-sample assignments, fitted models,
#' and misclassification error rates across simulation runs:
#' \item{assign_store}{array-to-sample assignments for each study design}
#' \item{model_store}{models for each combination of study designs, normalization methods, and classification methods}
#' \item{error_store}{misclassification error rates of the specified
#' validation method(s) for each combination of study designs,
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
#' biological.effect.nc <- biological.effect[!rownames(biological.effect) %in%
#'   ctrl.genes, ]
#' handling.effect.nc <- handling.effect[!rownames(handling.effect) %in% ctrl.genes, ]
#'
#' group.id <- substr(colnames(biological.effect.nc), 7, 7)
#'
#' # randomly split biological effect data into training and test set with
#' # equal number of endometrial and ovarian samples
#' biological.effect.train.ind <- colnames(biological.effect.nc)[c(sample(which(
#'   group.id == "E"), size = 64), sample(which(group.id == "V"), size = 64))]
#' biological.effect.test.ind <- colnames(biological.effect.nc)[!colnames(
#'   biological.effect.nc) %in% biological.effect.train.ind]
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
#' precision.results.flex <- precision.simulate.flex(seed = 1, N = 3,
#'   biological.effect.tr = biological.effect.nc.tr,
#'   biological.effect.te = biological.effect.nc.te,
#'   handling.effect.tr = handling.effect.nc.tr,
#'   handling.effect.te = handling.effect.nc.te,
#'   group.id.tr = substr(colnames(biological.effect.nc.tr), 7, 7),
#'   group.id.te = substr(colnames(biological.effect.nc.te), 7, 7),
#'   design.tr.list = c("PC-", "PC-", "STR", "STR"),
#'   design.te.list = c("PC-", "STR", "PC-", "STR"),
#'   norm.tr.list = c("NN", "QN", "NN", "QN"),
#'   norm.te.list = c("NN", "NN", "QN", "QN"),
#'   class.list = c("PAM", "LASSO"),
#'   batch.id.tr = list(1:40,
#'     41:64,
#'     (129:160) - 64,
#'     (161:192) - 64),
#'   batch.id.te = list(1:32, 33:64))
#' }
#'
"precision.simulate.flex" <- function(seed, N,
     biological.effect.tr, biological.effect.te,
     handling.effect.tr, handling.effect.te,
     group.id.tr, group.id.te,
     design.tr.list,
     design.te.list = NULL,
     norm.tr.list = c("NN", "QN"),
     norm.te.list = NULL,
     class.list = c("PAM", "LASSO"),
     valid.list = c("int", "ext.uh", "ext.sim.nuh"),
     batch.id.tr = NULL,
     batch.id.te = NULL,
     icombat = FALSE, isva = FALSE, iruv = FALSE,
     biological.effect.tr.ctrl = NULL,
     handling.effect.tr.ctrl = NULL,
     norm.tr.funcs = NULL,
     norm.te.funcs = NULL,
     class.funcs = NULL,
     pred.funcs = NULL){

  stopifnot(length(unique(design.tr.list)) > 1) # ensure create.storage() works
  stopifnot(design.tr.list %in% c("CC+", "CC-", "PC+", "PC-", "BLK", "STR"))
  stopifnot(valid.list %in% c("int", "ext.uh", "ext.sim.nuh"))
  if(is.null(design.te.list)){
   design.te.list <- design.tr.list
  }
  stopifnot(design.te.list %in% c("CC+", "CC-", "PC+", "PC-", "BLK", "STR"))
  stopifnot(length(design.tr.list) == length(design.te.list))

  if(is.null(norm.te.list)){
    norm.te.list <- norm.tr.list
  }
  stopifnot(length(norm.tr.list) == length(norm.te.list))

  # stopifnot(length(design.tr.list) == length(norm.tr.list))

  n.design <- length(design.tr.list)
  n.norm <- length(norm.tr.list) # length of norm.te.list should be the same
  n.class <- length(class.list)

  design.tr.list2 <- tolower(gsub("[-]", "2",
                                  gsub("[+]", "1", design.tr.list)))
  design.te.list2 <- tolower(gsub("[-]", "2",
                                  gsub("[+]", "1", design.te.list)))

  assign_store <- create.storage(
    design.list = union(design.tr.list2, design.te.list2),
    norm.list = "",
    class.list = "train",
    validating.sets = "assign")

  model_store <- create.storage(
    design.list = unique(design.tr.list2),
    norm.list = norm.tr.list,
    class.list = class.list,
    validating.sets = "model")

  error_store <- create.storage(
    design.list = paste(design.tr.list2,
                        design.te.list2,
                        sep = "-"),
    norm.list = paste0(norm.tr.list,
                      "-f", norm.te.list),
    class.list = class.list,
    validating.sets = valid.list)

  for(k in 1:N){ # each of the N simulation
    cat(k, "round seed used:", seed + k, "\n")

    for(dd in 1:length(design.tr.list)){

      dd.tr <- design.tr.list[dd]
      dd.tr2 <- tolower(gsub("[-]", "2",
                             gsub("[+]", "1", dd.tr)))

      if(!exists(paste0(dd.tr2, ".tr"))){
        # avoid repeatedly hybridizing the same training set

        cat("- 1. rehybridze sim.nuh training set -", dd.tr, "\n")

        if(dd.tr == "STR"){
          batch.id.tr.param <- paste0(", batch.id = , batch.id.tr")
        } else batch.id.tr.param <- ""

        eval(parse(text = paste0(dd.tr2, ".tr.ind <- ",
            "array.to.sample.assign.func(design = dd.tr,
            seed = seed + k, n = ncol(handling.effect.tr)",
            batch.id.tr.param, ")")))

        param <- ifelse(icombat, ", icombat = TRUE",
                        ifelse(isva, ", isva = TRUE",
                               ifelse(iruv, ", iruv = TRUE,
                                      biological.effect.ctrl = biological.effect.tr.ctrl,
                                      handling.effect.ctrl = handling.effect.tr.ctrl",
                                      "")))

        eval(parse(text = paste0(dd.tr2,
                   ".tr <- rehybridize(biological.effect = biological.effect.tr,
                   handling.effect = handling.effect.tr, group.id = group.id.tr,
                   array.to.sample.assign = ", dd.tr2,
                   ".tr.ind", param, ")")))

        eval(parse(text = paste0("assign_store['assign', ][['", dd.tr2,
                                 "']]$train.[[k]] <- list(",
                                 dd.tr2, ".tr.ind)")))

        if(isva){ # isva has special output from rehybridize()

          eval(parse(text = paste0("fsvaobj <- fsva(",
                                   dd.tr2, ".tr$trainData, ",
                                   dd.tr2, ".tr$trainMod, ",
                                   dd.tr2, ".tr$trainSV, biological.effect.te)")))

          eval(parse(text = paste0(dd.tr2,
                                   ".tr <- fsvaobj$db # train.sva")))
          eval(parse(text = "biological.effect.te2 <- fsvaobj$new # test.fsva"))

        } else biological.effect.te2 <- biological.effect.te

      }

      if(sum(valid.list %in% "ext.sim.nuh") == 1){
        # if simulated nonuniformly-handled data is selected as external validation

        dd.te <- design.te.list[dd]
        dd.te2 <- tolower(gsub("[-]", "2", gsub("[+]", "1", dd.te)))

        if(!exists(paste0(dd.te2, ".te"))){
          # avoid repeatedly hybridizing the same test set

          cat("- 1. rehybridze sim.nuh test set -", dd.te, "\n")

          if(dd.te == "STR"){
            batch.id.te.param <- paste0(", batch.id = , batch.id.te")
          } else batch.id.te.param <- ""

          eval(parse(text = paste0(dd.te2,
                     ".te.ind <- array.to.sample.assign.func(
                     design = dd.te, seed = seed + k,
                     n = ncol(handling.effect.te)", batch.id.te.param, ")")))

          eval(parse(text = paste0(dd.te2,
                     ".te <- rehybridize(biological.effect = biological.effect.te,
                     handling.effect = handling.effect.te,
                     group.id = group.id.te,
                     array.to.sample.assign = ", dd.te2, ".te.ind)")))

          eval(parse(text = paste0("assign_store['assign', ][['",
                                   dd.te2, "']]$test.[[k]] <- list(",
                                   dd.te2, ".te.ind)")))

          if(isva){ # isva has special output from rehybridize()

            # temporary paste0(dd.tr2, ".tr")
            eval(parse(text = paste0(dd.tr2,
                                     ".tr.temp <- rehybridize(biological.effect = biological.effect.tr,
                                     handling.effect = handling.effect.tr, group.id = group.id.tr,
                                     array.to.sample.assign = ", dd.tr2,
                                     ".tr.ind", param, ")")))

            eval(parse(text = paste0("fsvaobj <- fsva(",
                                     dd.tr2, ".tr.temp$trainData, ",
                                     dd.tr2, ".tr.temp$trainMod, ",
                                     dd.tr2, ".tr.temp$trainSV, ", dd.te2, ".te)")))

            ## only over-ride paste0(dd.tr2, ".tr") if "ext.uh" is not selected
            if(sum(valid.list %in% "ext.uh") == 1){
              eval(parse(text = paste0(dd.tr2,
                                       ".tr <- fsvaobj$db # train.sva")))
            }

            eval(parse(text = paste0(dd.te2, ".te <- fsvaobj$new # test.fsva")))

          }

        }

      }

      new.norm.tr.funcs <- switch.norm.funcs.flex(norm.list = norm.tr.list,
                                          norm.funcs = norm.tr.funcs)
      new.norm.te.funcs <- switch.norm.funcs.flex(norm.list = norm.te.list,
                                             norm.funcs = norm.te.funcs)

      for(xx in 1:length(norm.tr.list)){
      # for(norm.tr.met in norm.tr.list){
        norm.tr.met <- norm.tr.list[xx]
        norm.tr.met2 <- tolower(norm.tr.met)
        norm.tr.func <- new.norm.tr.funcs[xx] # take unique only

        norm.te.met <- norm.te.list[xx]
        norm.te.met2 <- tolower(norm.te.met)
        norm.te.func <- new.norm.te.funcs[xx] # take unique only

        if(!exists(paste0(dd.tr2, ".tr.", norm.tr.met2))){
          # avoid repeatedly processing the same training data

          cat("- 2. normalize + summarize sim.nuh.data training set -",
              dd.tr, "with", norm.tr.met, "\n")

          eval(parse(text = paste0("temp <- ", norm.tr.func,
                                   "(train = ", dd.tr2, ".tr)")))
          eval(parse(text = paste0(dd.tr2, ".tr.", norm.tr.met2,
                                   " <- temp$train.", norm.tr.met2)))
          eval(parse(text = paste0(dd.tr2, ".ref.", norm.tr.met2,
                                   " <- temp$ref.dis")))
          eval(parse(text = paste0(dd.tr2, ".tr.", norm.tr.met2,
                                   ".fin <- med.sum.pbset(", dd.tr2,
                                   ".tr.", norm.tr.met2, ")")))
          rm(temp)
        }

        if(sum(valid.list %in% "ext.uh") == 1){
          # if uniformly-handled data is selected as external validation

          if(!exists(paste0("uh.f", norm.te.met2, ".te.wrt.",
                            dd.tr2, ".", norm.tr.met2, ".tr"))){
            # avoid repeatedly processing the same uniformly-handled test data

            cat("- 2. normalize + summarize uh.data test set, with",
                norm.te.met, ", wrt", dd.tr, "\n")

            ## if norm.tr and norm.te are not the same, then rerun norm.tr.func based on norm.te.func... instead of creating ref.dis
            ## if norm.tr and norm.te are the same, then use ref.dis like current
            if(norm.tr.met2 == norm.te.met2){
              eval(parse(text = paste0("temp <- ", norm.te.func,
                                       "(test = biological.effect.te2",
                                       ", ref.dis = ", dd.tr2,
                                       ".ref.", norm.tr.met2, ")")))
            } else{
              eval(parse(text = paste0("temp <- ", norm.te.func,
                                       "(train = ", dd.tr2, ".tr.", norm.tr.met2,
                                       ", test = biological.effect.te2)")))
            }

            eval(parse(text = paste0("uh.f", norm.te.met2,
                              ".te.wrt.", dd.tr2, ".", norm.tr.met2, ".tr",
                              " <- temp$test.f", norm.te.met2)))
            eval(parse(text = paste0("uh.f", norm.te.met2,
                              ".te.wrt.", dd.tr2, ".", norm.tr.met2,
                              ".tr.fin <- med.sum.pbset(uh.f",
                              norm.te.met2, ".te.wrt.",
                              dd.tr2, ".", norm.tr.met2, ".tr)")))
            rm(temp)
          }
        }

        if(sum(valid.list %in% "ext.sim.nuh") == 1){
          # if simulated nonuniformly-handled data is selected as external validation

          if(!exists(paste0("sim.nuh.", dd.te2, ".f",
                            norm.te.met2, ".te.wrt.",
                            dd.tr2, ".", norm.tr.met2, ".tr"))){
              # avoid repeatedly processing the same simulated nonuniformly-handled test data

              dd.te <- design.te.list[dd]
              dd.te2 <- tolower(gsub("[-]", "2",
                                     gsub("[+]", "1", dd.te)))

              cat("- 2. normalize + summarize sim.nuh.data test set -",
                  dd.tr, "-", dd.te, "with", norm.te.met, "\n")

              if(norm.tr.met2 == norm.te.met2){
                eval(parse(text = paste0("temp <- ", norm.te.func,
                                  "(test = ", dd.te2, ".te",
                                  ", ref.dis = ", dd.tr2,
                                  ".ref.", norm.tr.met2, ")")))
              } else{
                eval(parse(text = paste0("temp <- ", norm.te.func,
                                         "(train = ", dd.tr2, ".tr.", norm.tr.met2,
                                         ", test = ", dd.te2, ".te)")))
              }

              eval(parse(text = paste0("sim.nuh.", dd.te2, ".f",
                                norm.te.met2, ".te.wrt.",
                                dd.tr2, ".", norm.tr.met2, ".tr",
                                " <- temp$test.f", norm.te.met2)))
              eval(parse(text = paste0("sim.nuh.", dd.te2, ".f",
                                norm.te.met2, ".te.wrt.",
                                dd.tr2, ".", norm.tr.met2,
                                ".tr.fin <- med.sum.pbset(",
                                "sim.nuh.", dd.te2, ".f",
                                norm.te.met2, ".te.wrt.",
                                dd.tr2, ".", norm.tr.met2, ".tr)")))
              rm(temp)
            }

        }

        new.class.funcs <- switch.classifier.funcs(
          class.list = class.list,
          class.funcs = class.funcs,
          pred.funcs = pred.funcs)

        for(cc in class.list){

          class.func <- new.class.funcs$build.funcs[class.list == cc]
          pred.func <- new.class.funcs$pred.funcs[class.list == cc]

          # build
          if(!exists(paste0(dd.tr2, ".", tolower(cc), ".int.", norm.tr.met2))){
            # avoid repeatedly building the same model

            cat(paste0("- 3. build ", cc, " ", dd.tr, " when ", norm.tr.met, "\n"))

            eval(parse(text = paste0(dd.tr2, ".",
                       tolower(cc), ".int.", norm.tr.met2, " <- ",
                       class.func, "(kfold = 5, X = ",
                       dd.tr2, ".tr.", norm.tr.met2,
                       ".fin, y = group.id.tr, seed = ",
                       seed + k, ")")))

            # store model
            eval(parse(text = paste0("model_store['model', ][['",
                                     dd.tr2, "']][['", cc, ".",
                                     norm.tr.met, "']][[k]] <- list(", dd.tr2, ".",
                                     tolower(cc), ".int.", norm.tr.met2, ")")))

            # store feature
            # can extract later
          }

          ## needs to repeatedly storing the same model and internal error

          # store internal error
          if(sum(valid.list %in% "int") == 1){
            for(cur.tr in which(design.tr.list == dd.tr)){
              eval(parse(text = paste0("error_store['int', ][[",
                           cur.tr, "]][['", cc, ".",
                           norm.tr.met, "-f", norm.te.met, "']][k] <- ", dd.tr2, ".",
                           tolower(cc), ".int.", norm.tr.met2, "$mc")))
            }
          }

          # predict
          if(sum(valid.list %in% "ext.uh") == 1){
            # if uniformly-handled data is selected as external validation

            if(!exists(paste0("uh.", dd.tr2, ".", norm.tr.met2, ".tr.",
                              tolower(cc), ".pred.", norm.te.met2))){
              # avoid repeatedly predicting with the same uniformly-handled data

              cat("- 4. predict based on", norm.te.met, "uh.data test data -",
                  dd.tr, "with", norm.tr.met, "and", cc, "\n")

              eval(parse(text = paste0("uh.", dd.tr2, ".", norm.tr.met2, ".tr.",
                           tolower(cc), ".pred.", norm.te.met2,
                           " <- ", pred.func, "(",
                           dd.tr2, ".", tolower(cc),
                           ".int.", norm.tr.met2, ", uh.f",
                           norm.te.met2, ".te.wrt.",
                           dd.tr2, ".", norm.tr.met2,
                           ".tr.fin, group.id.te)")))
            }

            # store external error
            for(cur.tr in which(design.tr.list == dd.tr)){
              eval(parse(text = paste0("error_store['ext.uh', ][[",
                                       cur.tr, "]][['", cc, ".",
                                       norm.tr.met, "-f", norm.te.met, "']][k] <- uh.",
                                       dd.tr2, ".", norm.tr.met2, ".tr.",
                                       tolower(cc), ".pred.", norm.te.met2, "$mc")))
            }

            print(eval(parse(text = paste0("uh.", dd.tr2, ".",
                                           norm.tr.met2, ".tr.",
                                           tolower(cc), ".pred.",
                                           norm.te.met2, "$mc"))))

          }

          if(sum(valid.list %in% "ext.sim.nuh") == 1){
            # if simulated nonuniformly-handled data is selected as external validation

            if(!exists(paste0("sim.nuh.", dd.te2, ".te.",
                              dd.tr2, ".", norm.tr.met2, ".tr.",
                              tolower(cc), ".pred.", norm.te.met2))){
              # avoid repeatedly predicting with the same simulated nonuniformly-handled data

              cat("- 4. predict based on", norm.te.met, "sim.nuh.data test data -",
                  dd.tr, "-", dd.te, "with", norm.tr.met, "and", cc, "\n")

              eval(parse(text = paste0("sim.nuh.", dd.te2, ".te.",
                     dd.tr2, ".", norm.tr.met2, ".tr.",
                     tolower(cc), ".pred.", norm.te.met2,
                     " <- ", pred.func, "(", dd.tr2, ".", tolower(cc),
                     ".int.", norm.tr.met2, ", sim.nuh.", dd.te2,
                     ".f", norm.te.met2, ".te.wrt.", dd.tr2, ".", norm.tr.met2,
                     ".tr.fin, group.id.te)")))
            }

            # store external error
            eval(parse(text = paste0("error_store['ext.sim.nuh', ][['",
                         dd.tr2, "-", dd.te2, "']][['", cc, ".",
                         norm.tr.met, "-f", norm.te.met, "']][k] <- sim.nuh.",
                         dd.te2, ".te.", dd.tr2, ".", norm.tr.met2, ".tr.",
                         tolower(cc), ".pred.", norm.te.met2, "$mc")))

            print(eval(parse(text = paste0("sim.nuh.",
                         dd.te2, ".te.", dd.tr2, ".", norm.tr.met2, ".tr.",
                         tolower(cc), ".pred.", norm.te.met2, "$mc"))))

          }

        } # class.list loop

      } # norm.list loop

    } # design.list loop

    # clean up objects before next simulation
    trycatch.func(rm(list = ls()[!ls() %in% c("seed", "k", "N",
                     "biological.effect.tr", "biological.effect.te", "biological.effect.te2",
                     "handling.effect.tr", "handling.effect.te",
                     "group.id.tr", "group.id.te",
                     "design.tr.list", "design.te.list",
                     "norm.tr.list", "norm.te.list",
                     "class.list", "valid.list",
                     "batch.id.tr", "batch.id.te",
                     "icombat", "isva", "iruv",
                     "biological.effect.tr.ctrl", "handling.effect.tr.ctrl",
                     "norm.tr.funcs", "norm.te.funcs",
                     "class.funcs", "pred.funcs",
                     "assign_store", "model_store", "error_store")]))

  } # simulation loop

  return(list(assign_store = assign_store,
              model_store = model_store,
              error_store = error_store))
}
