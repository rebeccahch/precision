## ----prestore.obj, echo = FALSE, error = FALSE, warning = FALSE, results = "hide"----
# load pre-stored R objects for the vignette.
load(file = "prestore.obj.Rdata")


## ----load.pkg, message = FALSE-------------------------------------------
if(!require(PRECISION)) install.packages("PRECISION")
library(PRECISION)


## ----load.data, eval = FALSE---------------------------------------------
#  data("uhdata.pl")
#  data("nuhdata.pl")
#  

## ----data.tab, results = "asis", echo = FALSE----------------------------
knitr::kable(uhdata.pl[1:15, 1:5], caption = "Uniformly-handled Data")

knitr::kable(nuhdata.pl[1:15, 1:5], caption = "nonuniformly-handled Data")


## ----per.unipbset.truncate, eval = FALSE---------------------------------
#  uhdata.pl.p5 <- per.unipbset.truncate(data = uhdata.pl,
#                                        num.per.unipbset = 5)

## ----ctrl.genes----------------------------------------------------------
# negative control probes
ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]

uhdata.pl.nc <- uhdata.pl[!rownames(uhdata.pl) %in% ctrl.genes, ] # nc for non-control probes

## ----est.eff, results = "hide", message = FALSE--------------------------
biological.effect <- estimate.biological.effect(uhdata = uhdata.pl)

handling.effect <- estimate.handling.effect(uhdata = uhdata.pl, 
                             nuhdata = nuhdata.pl)

biological.effect.nc <- biological.effect[!rownames(biological.effect) %in% ctrl.genes, ]
handling.effect.nc <- handling.effect[!rownames(handling.effect) %in% ctrl.genes, ]


## ----data.split, comment = ">", message = FALSE, eval = FALSE------------
#  set.seed(101)
#  group.id <- substr(colnames(biological.effect.nc), 7, 7)
#  
#  # randomly split biological effect data into training and test set with equal number of endometrial and ovarian samples
#  biological.effect.train.ind <- colnames(biological.effect.nc)[c(sample(which(group.id == "E"), size = 64), sample(which(group.id == "V"), size = 64))]
#  biological.effect.test.ind <- colnames(biological.effect.nc)[!colnames(biological.effect.nc) %in% biological.effect.train.ind]
#  
#  # non-randomly split handling effect data into training and test set
#  handling.effect.train.ind <- colnames(handling.effect.nc)[c(1:64, 129:192)]
#  handling.effect.test.ind <- colnames(handling.effect.nc)[65:128]
#  
#  group.id.list <- list("all" = group.id,
#                   "tr" = substr(biological.effect.train.ind, 7, 7),
#                   "te" = substr(biological.effect.test.ind, 7, 7))
#  array.to.sample.assign <- list("all" = c(rep(c("E", "V"), each = 64),
#                                      rep(c("V", "E"), each = 32)),
#                            "tr" = rep(c("E", "V"), each = 64),
#                            "te" = rep(c("V", "E"), each = 32))
#  

## ----study.design--------------------------------------------------------
set.seed(101)

# Complete Confounding
cc.ind <- confounding.design(seed = 1, num.array = 128, 
                               degree = "complete", rev.order = FALSE)

# Partial Confounding Reversed
pc.rev.ind <- confounding.design(seed = 1, num.array = 128, 
                               degree = "partial", rev.order = TRUE)

# Stratification
batch.id <- list(1:40, 41:64, (129:160) - 64, (161:192) - 64)
str.ind <- stratification.design(seed = 1, num.array = 128,
                                       batch.id = batch.id)

# Blocking
blk.ind <- blocking.design(seed = 1, num.array = 128)


## ----study.design.fig, fig.align = "hold", fig.height = 4, fig.width = 6, echo = FALSE----
par(mar = c(1, 1, 1, 1))
plot(1:128, pch = "",
     xlab = "", ylab = "", yaxt = "n", xaxt = "n", 
     ylim = c(0.25, 2.5))

points(1:128, rep(2, 128), 
     col = ifelse(cc.ind < 65, 2, 4), pch = "|")

points(1:128, rep(1.5, 128), 
     col = ifelse(pc.rev.ind < 65, 2, 4), pch = "|")

points(1:128, rep(1, 128), 
     col = ifelse(str.ind < 65, 2, 4), pch = "|")

points(1:128, rep(0.5, 128), 
     col = ifelse(blk.ind < 65, 2, 4), pch = "|")

text(x = 64.5, y = 2.25, labels = "Complete confounding")
text(x = 64.5, y = 1.75, labels = "Partial confounding reversed")
text(x = 64.5, y = 1.25, labels = "Stratification")
text(x = 64.5, y = 0.75, labels = "Blocking")


## ----rehybridize, eval = FALSE-------------------------------------------
#  assign.ind <- confounding.design(seed = 1, num.array = 192,
#                                 degree = "complete", rev.order = FALSE)
#  group.id <- substr(colnames(biological.effect.nc), 7, 7)
#  
#  # re-hybridize
#  sim.data.raw <- rehybridize(biological.effect = biological.effect.nc,
#                              handling.effect = handling.effect.nc,
#                              group.id = group.id,
#                              array.to.sample.assign = assign.ind)
#  
#  # re-hybridize + correct batch effects with SVA
#  sim.data.sva <- rehybridize(biological.effect = biological.effect.nc,
#                              handling.effect = handling.effect.nc,
#                              group.id = group.id,
#                              array.to.sample.assign = assign.ind,
#                              isva = TRUE)
#  
#  # re-hybridize + correct batch effects with RUV-4
#  biological.effect.ctrl <- biological.effect[rownames(biological.effect) %in% ctrl.genes, ]
#  handling.effect.ctrl <- handling.effect[rownames(handling.effect) %in% ctrl.genes, ]
#  
#  sim.data.ruv <- rehybridize(biological.effect = biological.effect.nc,
#                              handling.effect = handling.effect.nc,
#                              group.id = group.id,
#                              array.to.sample.assign = assign.ind,
#                              iruv = TRUE,
#                              biological.effect.ctrl = biological.effect.ctrl,
#                              handling.effect.ctrl = handling.effect.ctrl)
#  

## ----normalize, comment = ">", message = FALSE, eval = FALSE-------------
#  set.seed(101)
#  group.id <- substr(colnames(sim.data.raw), 7, 7)
#  
#  # randomly split data into training and test set with equal number of endometrial and ovarian samples
#  train.ind <- colnames(biological.effect)[c(sample(which(group.id == "E"), size = 64), sample(which(group.id == "V"), size = 64))]
#  
#  train.dat <- sim.data.raw[, train.ind]
#  test.dat <- sim.data.raw[, !colnames(sim.data.raw) %in%train.ind]
#  
#  
#  # median normalize (normalize training data only)
#  data.mn <- med.norm(train = train.dat)
#  
#  # median normalize (normalize training data + frozen normaliz test data)
#  data.mn <- med.norm(train = train.dat, test = test.dat)
#  
#  # quantile normalize
#  data.qn <- quant.norm(train = train.dat, test = test.dat)
#  
#  # varaince stabilizing normalize
#  data.vsn <- vs.norm(train = train.dat, test = test.dat)
#  

## ----med.sum.pbset, eval = FALSE-----------------------------------------
#  uhdata.psl <- med.sum.pbset(data = uhdata.pl,
#                              num.per.unipbset = 10)

## ----classification, eval = FALSE----------------------------------------
#  set.seed(101)
#  # randomly split biological effect data into training and test set with equal number of endometrial and ovarian samples
#  biological.effect.train.ind <- colnames(biological.effect.nc)[c(sample(which(group.id == "E"), size = 64), sample(which(group.id == "V"), size = 64))]
#  biological.effect.test.ind <- colnames(biological.effect.nc)[!colnames(biological.effect.nc) %in% biological.effect.train.ind]
#  
#  biological.effect.nc.tr <- biological.effect.nc[, biological.effect.train.ind]
#  biological.effect.nc.te <- biological.effect.nc[, biological.effect.test.ind]
#  
#  # build a PAM classifier
#  pam.int <- pam.intcv(X = biological.effect.nc.tr,
#                       y = substr(colnames(biological.effect.nc.tr), 7, 7),
#                       kfold = 5, seed = 1)
#  
#  # predict with the PAM classifier
#  pam.pred <- pam.predict(pam.intcv.model = pam.int,
#                          pred.obj = biological.effect.nc.te,
#                          pred.obj.group.id = substr(colnames(biological.effect.nc.te), 7, 7))
#  
#  pam.int$mc
#  pam.pred$mc
#  
#  # build a LASSO classifier
#  lasso.int <- lasso.intcv(X = biological.effect.nc.tr,
#                       y = substr(colnames(biological.effect.nc.tr), 7, 7),
#                       kfold = 5, seed = 1, alp = 1)
#  
#  # predict with the LASSO classifier
#  lasso.pred <- lasso.predict(lasso.intcv.model = lasso.int,
#                          pred.obj = biological.effect.nc.te,
#                          pred.obj.group.id = substr(colnames(biological.effect.nc.te), 7, 7))
#  
#  lasso.int$mc
#  lasso.pred$mc
#  

## ----uni.handled.simulate, results = "hide", eval = FALSE----------------
#  uni.handled.results <- uni.handled.simulate(seed = 1, N = 3,
#                                              biological.effect = biological.effect.nc,
#                                              norm.list = c("NN", "QN"),
#                                              class.list = c("PAM", "LASSO"))
#  

## ----uni.handled.tab, results = "asis"-----------------------------------
knitr::kable(data.frame(uni.handled.results$error_store[2, ]), 
             caption = "Analysis of uniformly-handled data")


## ----precision.simulate, results = "hide", eval = FALSE------------------
#  set.seed(101)
#  group.id <- substr(colnames(biological.effect.nc), 7, 7)
#  
#  # randomly split biological effect data into training and test set with equal number of endometrial and ovarian samples
#  biological.effect.train.ind <- colnames(biological.effect.nc)[c(sample(which(group.id == "E"), size = 64), sample(which(group.id == "V"), size = 64))]
#  biological.effect.test.ind <- colnames(biological.effect.nc)[!colnames(biological.effect.nc) %in% biological.effect.train.ind]
#  
#  biological.effect.train.test.split =
#    list("tr" = biological.effect.train.ind,
#         "te" = biological.effect.test.ind)
#  
#  # non-randomly split handling effect data into training and test set; technician effect as proxy
#  handling.effect.train.test.split =
#    list("tr" = c(1:64, 129:192),
#         "te" = 65:128)
#  
#  biological.effect.nc.tr <- biological.effect.nc[, biological.effect.train.ind]
#  biological.effect.nc.te <- biological.effect.nc[, biological.effect.test.ind]
#  handling.effect.nc.tr <- handling.effect.nc[, c(1:64, 129:192)]
#  handling.effect.nc.te <- handling.effect.nc[, 65:128]
#  
#  # Simulation without batch adjustment
#  precision.results <- precision.simulate(seed = 1, N = 3,
#                            biological.effect.tr = biological.effect.nc.tr,
#                            biological.effect.te = biological.effect.nc.te,
#                            handling.effect.tr = handling.effect.nc.tr,
#                            handling.effect.te = handling.effect.nc.te,
#                            group.id.tr = substr(colnames(biological.effect.nc.tr), 7, 7),
#                            group.id.te = substr(colnames(biological.effect.nc.te), 7, 7),
#                            design.list = c("PC-", "STR"),
#                            norm.list = c("NN", "QN"),
#                            class.list = c("PAM", "LASSO"),
#                            batch.id = list(1:40, 41:64, (129:160) - 64, (161:192) - 64))
#  
#  # Simulation with RUV-4 batch adjustment
#  biological.effect.tr.ctrl <- biological.effect.ctrl[, biological.effect.train.test.split$tr]
#  handling.effect.tr.ctrl <- handling.effect.ctrl[, handling.effect.train.test.split$tr]
#  
#  precision.ruv4.results <- precision.simulate(seed = 1, N = 3,
#                                biological.effect.tr = biological.effect.nc.tr,
#                                biological.effect.te = biological.effect.nc.te,
#                                handling.effect.tr = handling.effect.nc.tr,
#                                handling.effect.te = handling.effect.nc.te,
#                                group.id.tr = substr(colnames(biological.effect.nc.tr), 7, 7),
#                                group.id.te = substr(colnames(biological.effect.nc.te), 7, 7),
#                                design.list = c("PC-", "STR"),
#                                norm.list = c("NN", "QN"),
#                                class.list = c("PAM", "LASSO"),
#                                batch.id = list(1:40, 41:64, (129:160) - 64, (161:192) - 64),
#                                iruv = TRUE,
#                                biological.effect.tr.ctrl = biological.effect.tr.ctrl,
#                                handling.effect.tr.ctrl = handling.effect.tr.ctrl)
#  

## ----precision.tab, results = "asis"-------------------------------------
knitr::kable(data.frame(precision.results$error_store[2, "PC-"]), 
             aligh = "l",
             caption = "Analysis of simulated data in presence of confounding handling")

knitr::kable(data.frame(precision.ruv4.results$error_store[2, "PC-"]), 
             aligh = "l",
             caption = "Analysis of simulated data in presence of confounding handling, adjusting batch effects with RUV-4")


## ----norm.func.extend, eval = FALSE--------------------------------------
#  # an example of user-defined normalization method function:
#  "rand.norm" <- function(train, test = NULL){
#    stopifnot(nrow(train) == nrow(test))
#  
#    if(is.null(test)) {
#      test.frn <- NULL
#    } else{
#      train.rn <- apply(train, 2, sample)
#      dimnames(train.rn) <- dimnames(train)
#  
#      test.frn <- apply(test, 2, sample)
#      dimnames(test.frn) <- dimnames(test)
#    }
#  
#    return(list("train.rn" = train.rn,
#                "test.frn" = test.frn))
#  }
#  
#  uni.handled.results.rn <- uni.handled.simulate(seed = 1, N = 3,
#                                         biological.effect = biological.effect.nc,
#                                         norm.list = c("NN", "QN", "RN"),
#                                         class.list = c("PAM", "LASSO"),
#                                         norm.funcs = "rand.norm")
#  

## ----uni.handled.tab.rn, results = "asis"--------------------------------
knitr::kable(data.frame(uni.handled.results.rn$error_store[2, ]), 
             caption = "Analysis of uniformly-handled data (with RN)")


## ----class.func.extend, eval = FALSE-------------------------------------
#  # example of user-defined classification method functions:
#  "build.ridge" <- function(kfold = 5, X, y, seed, alp = 0){
#    ptm <- proc.time()
#    set.seed(seed)
#  
#    cv.fit <- glmnet::cv.glmnet(x = data.matrix(t(X)), y = factor(y),
#                                family = "binomial", type.measure = "class", alpha = alp, nfold = kfold)
#    mc <- cv.fit$cvm[which(cv.fit$lambda == cv.fit$lambda.1se)]
#    #best.lambda <- cv.fit$lambda.1se # can be extracted from cv.fit
#    coefs <- coef(cv.fit, s = "lambda.1se")
#    time <- proc.time() - ptm
#    return(list(mc=mc, time=time, model=cv.fit, cfs=coefs))
#  }
#  
#  # Ridge has the same prediction function as LASSO
#  
#  uni.handled.results.ridge <- uni.handled.simulate(seed = 1, N = 3,
#                                     biological.effect = biological.effect.nc,
#                                     norm.list = c("NN", "RN"),
#                                     class.list = c("LASSO", "Ridge"),
#                                     norm.funcs = "rand.norm",
#                                     class.funcs = "build.ridge",
#                                     pred.funcs = "lasso.predict")
#  

## ----uni.handled.ridge.tab, results = "asis"-----------------------------
knitr::kable(data.frame(uni.handled.results.ridge$error_store[2, ]), 
             caption = "Analysis of uniformly-handled data (with Ridge)")


## ----dea, eval = FALSE---------------------------------------------------
#  uhdata.psl.nc <- uhdata.psl[!rownames(uhdata.psl) %in% ctrl.genes, ] # nc for non-control probes
#  
#  group.id <- substr(colnames(uhdata.psl.nc), 7, 7)
#  group.id.level <- levels(as.factor(group.id))
#  limma.fit.uhdata<- limma.pbset(data = uhdata.psl.nc,
#                                 group.id = group.id,
#                                 group.id.level = group.id.level)
#  table(limma.fit.uhdata$P.Value < 0.01, dnn = "DE genes")
#  

## ----dea.tab, results = "asis"-------------------------------------------
tab <- data.frame(table(limma.fit.uhdata$P.Value < 0.01))
colnames(tab) <- c("DEA", "Count")
knitr::kable(tab, rownames = NULL)


## ----reduce.signal, eval = FALSE-----------------------------------------
#  # reduced signal by half
#  group.id <- substr(colnames(biological.effect.nc), 7, 7)
#  redhalf.biological.effect.nc <- reduce.signal(biological.effect = biological.effect.nc,
#                              group.id = group.id,
#                              group.id.level = c("E", "V"),
#                              reduce.multiplier = 1/2)
#  
#  
#  # extract differential expressed genes
#  biological.effect.nc.psl <- med.sum.pbset(biological.effect.nc)
#  s.e.limma.fit <- limma.pbset(data = biological.effect.nc.psl,
#                           group.id = group.id,
#                           group.id.level = c("E", "V"))
#  de.ind <- s.e.limma.fit$P.Value < 0.01
#  
#  de.gene.name <- rownames(redhalf.biological.effect.nc)[which(de.ind)][2]
#  nonde.gene.name <- rownames(redhalf.biological.effect.nc)[-which(de.ind)][2]
#  
#  biological.effect.nc.de <- biological.effect.nc[de.gene.name, ]
#  redhalf.biological.effect.nc.de <- redhalf.biological.effect.nc[de.gene.name, ]
#  
#  biological.effect.nc.nonde <- biological.effect.nc[nonde.gene.name, ]
#  redhalf.biological.effect.nc.nonde <- redhalf.biological.effect.nc[nonde.gene.name, ]
#  

## ----reduce.signal.fig, echo = FALSE, fig.align = "hold", fig.width = 4, fig.height = 4----
group.id <- substr(colnames(biological.effect.nc), 7, 7)

# plot
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))
# for DE genes, the difference has been shrunken in half...
boxplot(biological.effect.nc.de ~ group.id == "E", main = "DE")
boxplot(redhalf.biological.effect.nc.de ~ group.id == "E", main = "DE")
# for non-DE genes, the difference has not been changed
boxplot(biological.effect.nc.nonde ~ group.id == "E", main = "Non-DE")
boxplot(redhalf.biological.effect.nc.nonde ~ group.id == "E", main = "Non-DE")


## ----amplify.handling.effect, fig.align = "hold", fig.width = 7, fig.height = 7----
handling.effect.nc.tr <- handling.effect.nc[, c(1:64, 129:192)]

# shift
handling.effect.nc.tr.shift <- amplify.handling.effect(handling.effect = handling.effect.nc.tr,
       amplify.array.id = colnames(handling.effect.nc.tr)[1:64],
       amplify.level = 2, type = "shift")

# scale 1
handling.effect.nc.tr.scale1 <- amplify.handling.effect(handling.effect = handling.effect.nc.tr,
       amplify.array.id = colnames(handling.effect.nc.tr)[1:64],
       amplify.level = 2, type = "scale1")

# scale 2
amplify.array.id <- list(1:40, 41:64, (129:160) - 64, (161:192) - 64)
for(i in 1:length(amplify.array.id)) 
  amplify.array.id[[i]] <- colnames(handling.effect.nc.tr)[amplify.array.id[[i]]]
amplify.level <- c(1.2, 1.3, 1/3, 2/3)

handling.effect.nc.tr.scale2 <- amplify.handling.effect(handling.effect = handling.effect.nc.tr,
       amplify.array.id = amplify.array.id,
       amplify.level = amplify.level,
       type = "scale2")

par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))
rng <- range(handling.effect.nc.tr, handling.effect.nc.tr.shift, 
             handling.effect.nc.tr.scale1, handling.effect.nc.tr.scale2)
boxplot(handling.effect.nc.tr, main = "Original",
        ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
boxplot(handling.effect.nc.tr.shift, main = "Shift",
        ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
boxplot(handling.effect.nc.tr.scale1, main = "Scaling 1",
        ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
boxplot(handling.effect.nc.tr.scale2, main = "Scaling 2",
        ylim = rng, pch = 20, cex = 0.2, xaxt = "n")


## ----calc.confounding.level, comment = ">", message = FALSE, eval = FALSE----
#  # calculate confounding level
#  nbe.genes <- ifelse(gene.cat == -1, TRUE, FALSE)
#  calc.confounding.level(data = biological.effect.nc[, biological.effect.train.ind],
#                         group.id = substr(biological.effect.train.ind, 7, 7),
#                         nbe.genes = nbe.genes)
#  

