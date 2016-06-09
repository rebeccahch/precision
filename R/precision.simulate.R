#' Classification analysis of simulation study
#'
#' Performs simulation study in Qin et al. (see reference).
#'
#' @references http://clincancerres.aacrjournals.org/content/20/13/3371.long
#' @param myseed specifies seed for random assignment using set.seed().
#' @param N number of simulation runs.
#' @param smp.eff.tr sample effect training data, rows as probes, columns as samples.
#' @param smp.eff.te sample effect test data, rows as probes, columns as samples; must have same number of probes and probe names as sample effect training data.
#' @param ary.eff.tr array effect training data, rows as probes, columns as samples; must have same dimensions and same probe name as sample effect training data.
#' @param ary.eff.te array effect test data, rows as probes, columns as samples; must have same dimensions and same probe name as array effect test data.
#' @param group.id.tr sample group ID of the sample effect training data; has to be "E" and "V".
#' @param group.id.te sample group ID of the sample effect test data; has to be "E" and "V".
#' @param design.list a list of strings for study designs compared in the simulation study;
#' built-in designs are "CC+", "CC-", "PC+", "PC-", "BLK", and "STR" for "Complete Confounding 1", "Complete Confounding 2", "Partial Confounding 1", "Partial Confounding 2", "Blocking", "Stratification" in Qin et al.
#' @param norm.list a list of strings for normalization methods compared in the simulation study;
#' build-in available normalization methods are "NN", "QN", "MN", "VSN" for "No Normalization", "Quantile Normalization", "Median Normalization", "Variance Stablizing Normalization";
#' user can provide a list of normalization methods given the functions are supplied (also see norm.funcs).
#' @param class.list a list of strings for classification methods compared in the simulation study;
#' built-in classification methods are "PAM" and "LASSO" for "prediction analysis for microarrays" and "least absolute shrinkage and selection operator";
#' user can provide a list of classification methods given the correponding model-building and predicting functions are supplied (also see class.funcs and pred.funcs).
#' @param batch.id a list of batch id by the number of batches when stratification study design is specified; default is NULL.
#' @param icombat indicator for combat adjustment; default is not to adjust (icombat = FALSE).
#' @param isva indicator for sva adjustment; default is not to adjust (isva = FALSE).
#' @param iruv indicator for RUV-4 adjustment; default is not to adjust (iruv = FALSE).
#' @param smp.eff.tr.ctrl negative-control gene sample effect data if iruv = TRUE, rows as probes, columns as samples;
#' must have same number of probes and probe names as non-control gene sample effect training data.
#' @param ary.eff.tr.ctrl negative-control gene array effect data if iruv = TRUE, rows as probes, columns as samples;
#' must have same number of probes and probe names as non-control gene array effect training data.
#' @param norm.funcs a list of strings for names of user-defined normalization method functions, in the order of norm.list excluding any built-in normalization methods.
#' @param class.funcs a list of strings for names of user-defined classification model-building functions, in the order of class.list excluding any built-in classification methods.
#' @param pred.funcs a list of strings for names of user-defined classification predicting functions, in the order of class.list excluding any built-in classification methods.
#' @return simulated results with list of models built and internal and external misclassification error stored, also a list of assignment stored
#' @keywords simulation
#' @export
#' @examples
#' \dontrun{
#' set.seed(101)
#' smp.eff <- estimate.smp.eff(r.data = r.data.pl)
#' ary.eff <- estimate.ary.eff(r.data = r.data.pl,
#'                              non.r.data = non.r.data.pl)
#'
#' ctrl.genes <- unique(rownames(r.data.pl))[grep("NC", unique(rownames(r.data.pl)))]
#'
#' smp.eff.nc <- smp.eff[!rownames(smp.eff) %in% ctrl.genes, ]
#' ary.eff.nc <- ary.eff[!rownames(ary.eff) %in% ctrl.genes, ]
#'
#' group.id <- substr(colnames(smp.eff.nc), 7, 7)
#'
#' # randomly split sample effect data into training and test set with
#' # equal number of endometrial and ovarian samples
#' smp.eff.train.ind <- colnames(smp.eff.nc)[c(sample(which(group.id == "E"), size = 64),
#'                                           sample(which(group.id == "V"), size = 64))]
#' smp.eff.test.ind <- colnames(smp.eff.nc)[!colnames(smp.eff.nc) %in% smp.eff.train.ind]
#' smp.eff.train.test.split =
#'   list("tr" = smp.eff.train.ind,
#'        "te" = smp.eff.test.ind)
#'
#' # non-randomly split array effect data into training and test set
#' ary.eff.train.test.split =
#'   list("tr" = c(1:64, 129:192),
#'        "te" = 65:128)
#'
#' smp.eff.nc.tr <- smp.eff.nc[, smp.eff.train.ind]
#' smp.eff.nc.te <- smp.eff.nc[, smp.eff.test.ind]
#' ary.eff.nc.tr <- ary.eff.nc[, c(1:64, 129:192)]
#' ary.eff.nc.te <- ary.eff.nc[, 65:128]
#'
#' # Simulation without batch adjustment
#' precision.results <- precision.simulate(myseed = 1, N = 3,
#'                                         smp.eff.tr = smp.eff.nc.tr,
#'                                         smp.eff.te = smp.eff.nc.te,
#'                                         ary.eff.tr = ary.eff.nc.tr,
#'                                         ary.eff.te = ary.eff.nc.te,
#'                                         group.id.tr = substr(colnames(smp.eff.nc.tr), 7, 7),
#'                                         group.id.te = substr(colnames(smp.eff.nc.te), 7, 7),
#'                                         design.list = c("PC-", "STR"),
#'                                         norm.list = c("NN", "QN"),
#'                                         class.list = c("PAM", "LASSO"),
#'                                         batch.id = list(1:40,
#'                                                         41:64,
#'                                                         (129:160) - 64,
#'                                                         (161:192) - 64))
#'
#' # Simulation with RUV-4 batch adjustment
#' smp.eff.ctrl <- smp.eff[rownames(smp.eff) %in% ctrl.genes, ]
#' ary.eff.ctrl <- ary.eff[rownames(ary.eff) %in% ctrl.genes, ]
#'
#' smp.eff.tr.ctrl <- smp.eff.ctrl[, smp.eff.train.test.split$tr]
#' ary.eff.tr.ctrl <- ary.eff.ctrl[, ary.eff.train.test.split$tr]
#'
#' precision.ruv4.results <- precision.simulate(myseed = 1, N = 3,
#'                                              smp.eff.tr = smp.eff.nc.tr,
#'                                              smp.eff.te = smp.eff.nc.te,
#'                                              ary.eff.tr = ary.eff.nc.tr,
#'                                              ary.eff.te = ary.eff.nc.te,
#'                                              group.id.tr = substr(colnames(smp.eff.nc.tr), 7, 7),
#'                                              group.id.te = substr(colnames(smp.eff.nc.te), 7, 7),
#'                                              design.list = c("PC-", "STR"),
#'                                              norm.list = c("NN", "QN"),
#'                                              class.list = c("PAM", "LASSO"),
#'                                              batch.id = list(1:40,
#'                                                              41:64,
#'                                                              (129:160) - 64,
#'                                                              (161:192) - 64),
#'                                              iruv = TRUE,
#'                                              smp.eff.tr.ctrl = smp.eff.tr.ctrl,
#'                                              ary.eff.tr.ctrl = ary.eff.tr.ctrl)
#' }

"precision.simulate" <- function(myseed, N,
                                 smp.eff.tr, smp.eff.te,
                                 ary.eff.tr, ary.eff.te,
                                 group.id.tr, group.id.te,
                                 design.list = c("CC+", "CC-", "PC+", "PC-"),
                                 norm.list = c("NN", "QN"),
                                 class.list = c("PAM", "LASSO"),
                                 batch.id = NULL,
                                 icombat = FALSE, isva = FALSE, iruv = FALSE,
                                 smp.eff.tr.ctrl = NULL,
                                 ary.eff.tr.ctrl = NULL,
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
    cat(k, "round seed used:", myseed + k, "\n")
    #cat("- setup simulated data \n")
    for(dd in design.list){
      if(dd == "CC+") cc1.ind <- confounding.design(seed = myseed + k, num.smp = ncol(ary.eff.tr), degree = "complete", rev.order = FALSE)
      if(dd == "CC-") cc2.ind <- confounding.design(seed = myseed + k, num.smp = ncol(ary.eff.tr), degree = "complete", rev.order = TRUE)
      if(dd == "PC+") pc1.ind <- confounding.design(seed = myseed + k, num.smp = ncol(ary.eff.tr), degree = "partial", rev.order = FALSE)
      if(dd == "PC-") pc2.ind <- confounding.design(seed = myseed + k, num.smp = ncol(ary.eff.tr), degree = "partial", rev.order = TRUE)
      if(dd == "BLK") blk.ind <- blocking.design(seed = myseed + k, num.smp = ncol(ary.eff.tr))
      if(dd == "STR") str.ind <- stratification.design(seed = myseed + k, num.smp = ncol(ary.eff.tr), batch.id = batch.id)

      dd2 <- tolower(gsub("[-]", "2", gsub("[+]", "1", dd)))

      param <- ifelse(icombat, ", icombat = TRUE",
                      ifelse(isva, ", isva = TRUE",
                             ifelse(iruv, ", iruv = TRUE,
                                    smp.eff.ctrl = smp.eff.tr.ctrl,
                                    ary.eff.ctrl = ary.eff.tr.ctrl", "")))

      eval(parse(text = paste0(dd2, ".tr <- rehybridize(smp.eff = smp.eff.tr,
                               ary.eff = ary.eff.tr, group.id = group.id.tr,
                               ary.to.smp.assign = ", dd2, ".ind", param, ")")))

      eval(parse(text = paste0("assign_store['assign', ][['", dd,
                               "']]$train.[[k]] <- list(", dd2, ".ind)")))

      if(isva){ # isva has special output from rehybridize()

        eval(parse(text = paste0("fsvaobj <- fsva(", dd2, ".tr$trainData, ",
                                 dd2, ".tr$trainMod, ",
                                 dd2, ".tr$trainSV, smp.eff.te)")))

        eval(parse(text = paste0(dd2, ".tr <- fsvaobj$db # train.sva")))
        eval(parse(text = "smp.eff.te2 <- fsvaobj$new # test.fsva"))
      } else {
        smp.eff.te2 <- smp.eff.te
      }

      #cat("- preprocess data \n")
      new.norm.funcs <- switch.norm.funcs(norm.list = norm.list, norm.funcs = norm.funcs)
      for(norm.met in norm.list){
        norm.met2 <- tolower(norm.met)
        norm.func <- new.norm.funcs[norm.list == norm.met]

        if(norm.met != "NN"){
          # normalize
          eval(parse(text = paste0("temp <- ", norm.func,"(", dd2, ".tr, smp.eff.te2)")))
          eval(parse(text = paste0(dd2, ".tr.", norm.met2, " <- temp$train.", norm.met2)))
          eval(parse(text = paste0("gs.", dd2, ".f", norm.met2, " <- temp$test.f", norm.met2)))

          # summarize
          eval(parse(text = paste0(dd2, ".tr.", norm.met2, ".fin <- med.sum.pbset(", dd2, ".tr.", norm.met2, ")")))
          eval(parse(text = paste0("gs.", dd2, ".f", norm.met2, ".fin <- med.sum.pbset(gs.", dd2, ".f", norm.met2, ")")))

        } else{
          eval(parse(text = paste0(dd2, ".tr.nn.fin <- med.sum.pbset(", dd2, ".tr)")))
          #eval(parse(text = paste0("gs.", dd2, ".fnn.fin <- med.sum.pbset(smp.eff.te2)")))
          eval(parse(text = paste0("temp <- quant.norm(", dd2, ".tr, smp.eff.te2)")))
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
                                   seed = ", myseed + k, ")")))
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
