## ----prestore.obj, echo = FALSE, error = FALSE, warning = FALSE, results = "hide"----
# load pre-stored R objects for the vignette.
load(file = "prestore.obj.Rdata")


## ----load.pkg, message = FALSE-------------------------------------------
if(!require(PRECISION)) install.packages("PRECISION")
library(PRECISION)


## ----load.data, eval = FALSE---------------------------------------------
#  data("r.data.pl")
#  data("non.r.data.pl")
#  

## ----data.tab, results = "asis", echo = FALSE----------------------------
knitr::kable(r.data.pl[1:15, 1:5], caption = "Uniformly-handled Data")

knitr::kable(non.r.data.pl[1:15, 1:5], caption = "Non-uniformly-handled Data")


## ----per.unipbset.truncate, eval = FALSE---------------------------------
#  r.data.pl.p5 <- per.unipbset.truncate(data = r.data.pl,
#                                        num.per.unipbset = 5)

## ----ctrl.genes----------------------------------------------------------
# negatively biological control probes
ctrl.genes <- unique(rownames(r.data.pl))[grep("NC", unique(rownames(r.data.pl)))]

r.data.pl.nc <- r.data.pl[!rownames(r.data.pl) %in% ctrl.genes, ] # nc for non-control probes

## ----est.eff, results = "hide", message = FALSE--------------------------
smp.eff <- estimate.smp.eff(r.data = r.data.pl)

ary.eff <- estimate.ary.eff(r.data = r.data.pl, 
                             non.r.data = non.r.data.pl)

smp.eff.nc <- smp.eff[!rownames(smp.eff) %in% ctrl.genes, ]
ary.eff.nc <- ary.eff[!rownames(ary.eff) %in% ctrl.genes, ]


## ----data.split, comment = ">", message = FALSE, eval = FALSE------------
#  set.seed(101)
#  group.id <- substr(colnames(smp.eff.nc), 7, 7)
#  
#  # randomly split sample effect data into training and test set with equal number of endometrial and ovarian samples
#  smp.eff.train.ind <- colnames(smp.eff.nc)[c(sample(which(group.id == "E"), size = 64), sample(which(group.id == "V"), size = 64))]
#  smp.eff.test.ind <- colnames(smp.eff.nc)[!colnames(smp.eff.nc) %in% smp.eff.train.ind]
#  
#  # non-randomly split array effect data into training and test set; technician effect as proxy
#  ary.eff.train.ind <- colnames(ary.eff.nc)[c(1:64, 129:192)]
#  ary.eff.test.ind <- colnames(ary.eff.nc)[65:128]
#  
#  group.id.list <- list("all" = group.id,
#                   "tr" = substr(smp.eff.train.ind, 7, 7),
#                   "te" = substr(smp.eff.test.ind, 7, 7))
#  ary.to.smp.assign <- list("all" = c(rep(c("E", "V"), each = 64),
#                                      rep(c("V", "E"), each = 32)),
#                            "tr" = rep(c("E", "V"), each = 64),
#                            "te" = rep(c("V", "E"), each = 32))
#  

## ----study.design--------------------------------------------------------
set.seed(101)

# Complete Confounding
cc.ind <- confounding.design(seed = 1, num.smp = 128, 
                               degree = "complete", rev.order = FALSE)

# Partial Confounding Reversed
pc.rev.ind <- confounding.design(seed = 1, num.smp = 128, 
                               degree = "partial", rev.order = TRUE)

# Stratification
batch.id <- list(1:40, 41:64, (129:160) - 64, (161:192) - 64)
str.ind <- stratification.design(seed = 1, num.smp = 128,
                                       batch.id = batch.id)

# Blocking
blk.ind <- blocking.design(seed = 1, num.smp = 128)


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
#  assign.ind <- confounding.design(seed = 1, num.smp = 192,
#                                 degree = "complete", rev.order = FALSE)
#  group.id <- substr(colnames(smp.eff.nc), 7, 7)
#  
#  # rehybridize
#  sim.data.raw <- rehybridize(smp.eff = smp.eff.nc,
#                              ary.eff = ary.eff.nc,
#                              group.id = group.id,
#                              ary.to.smp.assign = assign.ind)
#  
#  # rehybridize + adjust batch effects with SVA
#  sim.data.sva <- rehybridize(smp.eff = smp.eff.nc,
#                              ary.eff = ary.eff.nc,
#                              group.id = group.id,
#                              ary.to.smp.assign = assign.ind,
#                              isva = TRUE)
#  
#  # rehybridize + adjust batch effects with RUV-4
#  smp.eff.ctrl <- smp.eff[rownames(smp.eff) %in% ctrl.genes, ]
#  ary.eff.ctrl <- ary.eff[rownames(ary.eff) %in% ctrl.genes, ]
#  
#  sim.data.ruv <- rehybridize(smp.eff = smp.eff.nc,
#                              ary.eff = ary.eff.nc,
#                              group.id = group.id,
#                              ary.to.smp.assign = assign.ind,
#                              iruv = TRUE,
#                              smp.eff.ctrl = smp.eff.ctrl,
#                              ary.eff.ctrl = ary.eff.ctrl)
#  

## ----normalize, comment = ">", message = FALSE, eval = FALSE-------------
#  set.seed(101)
#  group.id <- substr(colnames(sim.data.raw), 7, 7)
#  
#  # randomly split data into training and test set with equal number of endometrial and ovarian samples
#  train.ind <- colnames(smp.eff)[c(sample(which(group.id == "E"), size = 64), sample(which(group.id == "V"), size = 64))]
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
#  r.data.psl <- med.sum.pbset(data = r.data.pl,
#                              num.per.unipbset = 10)

## ----classification, eval = FALSE----------------------------------------
#  set.seed(101)
#  # randomly split sample effect data into training and test set with equal number of endometrial and ovarian samples
#  smp.eff.train.ind <- colnames(smp.eff.nc)[c(sample(which(group.id == "E"), size = 64), sample(which(group.id == "V"), size = 64))]
#  smp.eff.test.ind <- colnames(smp.eff.nc)[!colnames(smp.eff.nc) %in% smp.eff.train.ind]
#  
#  smp.eff.nc.tr <- smp.eff.nc[, smp.eff.train.ind]
#  smp.eff.nc.te <- smp.eff.nc[, smp.eff.test.ind]
#  
#  # build a PAM classifier
#  pam.int <- pam.intcv(X = smp.eff.nc.tr,
#                       y = substr(colnames(smp.eff.nc.tr), 7, 7),
#                       kfold = 5, seed = 1)
#  
#  # predict with the PAM classifier
#  pam.pred <- pam.predict(pam.intcv.model = pam.int,
#                          pred.obj = smp.eff.nc.te,
#                          pred.obj.group.id = substr(colnames(smp.eff.nc.te), 7, 7))
#  
#  pam.int$mc
#  pam.pred$mc
#  
#  # build a LASSO classifier
#  lasso.int <- lasso.intcv(X = smp.eff.nc.tr,
#                       y = substr(colnames(smp.eff.nc.tr), 7, 7),
#                       kfold = 5, seed = 1, alp = 1)
#  
#  # predict with the LASSO classifier
#  lasso.pred <- lasso.predict(lasso.intcv.model = lasso.int,
#                          pred.obj = smp.eff.nc.te,
#                          pred.obj.group.id = substr(colnames(smp.eff.nc.te), 7, 7))
#  
#  lasso.int$mc
#  lasso.pred$mc
#  

## ----uni.handled.simulate, results = "hide", eval = FALSE----------------
#  uni.handled.results <- uni.handled.simulate(seed = 1, N = 3,
#                                              smp.eff = smp.eff.nc,
#                                              norm.list = c("NN", "QN"),
#                                              class.list = c("PAM", "LASSO"))
#  

## ----uni.handled.tab, results = "asis"-----------------------------------
knitr::kable(data.frame(uni.handled.results$error_store[2, ]), 
             caption = "Analysis of uniformly-handled data")


## ----precision.simulate, results = "hide", eval = FALSE------------------
#  set.seed(101)
#  group.id <- substr(colnames(smp.eff.nc), 7, 7)
#  
#  # randomly split sample effect data into training and test set with equal number of endometrial and ovarian samples
#  smp.eff.train.ind <- colnames(smp.eff.nc)[c(sample(which(group.id == "E"), size = 64), sample(which(group.id == "V"), size = 64))]
#  smp.eff.test.ind <- colnames(smp.eff.nc)[!colnames(smp.eff.nc) %in% smp.eff.train.ind]
#  
#  smp.eff.train.test.split =
#    list("tr" = smp.eff.train.ind,
#         "te" = smp.eff.test.ind)
#  
#  # non-randomly split array effect data into training and test set; technician effect as proxy
#  ary.eff.train.test.split =
#    list("tr" = c(1:64, 129:192),
#         "te" = 65:128)
#  
#  smp.eff.nc.tr <- smp.eff.nc[, smp.eff.train.ind]
#  smp.eff.nc.te <- smp.eff.nc[, smp.eff.test.ind]
#  ary.eff.nc.tr <- ary.eff.nc[, c(1:64, 129:192)]
#  ary.eff.nc.te <- ary.eff.nc[, 65:128]
#  
#  # Simulation without batch adjustment
#  precision.results <- precision.simulate(seed = 1, N = 3,
#                            smp.eff.tr = smp.eff.nc.tr,
#                            smp.eff.te = smp.eff.nc.te,
#                            ary.eff.tr = ary.eff.nc.tr,
#                            ary.eff.te = ary.eff.nc.te,
#                            group.id.tr = substr(colnames(smp.eff.nc.tr), 7, 7),
#                            group.id.te = substr(colnames(smp.eff.nc.te), 7, 7),
#                            design.list = c("PC-", "STR"),
#                            norm.list = c("NN", "QN"),
#                            class.list = c("PAM", "LASSO"),
#                            batch.id = list(1:40, 41:64, (129:160) - 64, (161:192) - 64))
#  
#  # Simulation with RUV-4 batch adjustment
#  smp.eff.tr.ctrl <- smp.eff.ctrl[, smp.eff.train.test.split$tr]
#  ary.eff.tr.ctrl <- ary.eff.ctrl[, ary.eff.train.test.split$tr]
#  
#  precision.ruv4.results <- precision.simulate(seed = 1, N = 3,
#                                smp.eff.tr = smp.eff.nc.tr,
#                                smp.eff.te = smp.eff.nc.te,
#                                ary.eff.tr = ary.eff.nc.tr,
#                                ary.eff.te = ary.eff.nc.te,
#                                group.id.tr = substr(colnames(smp.eff.nc.tr), 7, 7),
#                                group.id.te = substr(colnames(smp.eff.nc.te), 7, 7),
#                                design.list = c("PC-", "STR"),
#                                norm.list = c("NN", "QN"),
#                                class.list = c("PAM", "LASSO"),
#                                batch.id = list(1:40, 41:64, (129:160) - 64, (161:192) - 64),
#                                iruv = TRUE,
#                                smp.eff.tr.ctrl = smp.eff.tr.ctrl,
#                                ary.eff.tr.ctrl = ary.eff.tr.ctrl)
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
#                                         smp.eff = smp.eff.nc,
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
#                                     smp.eff = smp.eff.nc,
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
#  r.data.psl.nc <- r.data.psl[!rownames(r.data.psl) %in% ctrl.genes, ] # nc for non-control probes
#  
#  group.id <- substr(colnames(r.data.psl.nc), 7, 7)
#  group.id.level <- levels(as.factor(group.id))
#  limma.fit.r.data<- limma.pbset(data = r.data.psl.nc,
#                                 group.id = group.id,
#                                 group.id.level = group.id.level)
#  table(limma.fit.r.data$P.Value < 0.01, dnn = "DE genes")
#  

## ----dea.tab, results = "asis"-------------------------------------------
tab <- data.frame(table(limma.fit.r.data$P.Value < 0.01))
colnames(tab) <- c("DEA", "Count")
knitr::kable(tab, rownames = NULL)


## ----classify.gene.type, comment = ">", message = FALSE, eval = FALSE----
#  # classify gene type
#  gene.cat <- classify.gene.type(smp.eff = smp.eff.nc,
#                                 ary.eff = ary.eff.nc,
#                                 smp.eff.train.ind = smp.eff.train.ind,
#                                 ary.eff.train.ind = ary.eff.train.ind,
#                                 group.id = group.id.list,
#                                 ary.to.smp.assign = ary.to.smp.assign)
#  

## ----calc.confounding.level, comment = ">", message = FALSE, eval = FALSE----
#  # calculate confounding level
#  nbe.genes <- ifelse(gene.cat == -1, TRUE, FALSE)
#  calc.confounding.level(data = smp.eff.nc[, smp.eff.train.ind],
#                         group.id = substr(smp.eff.train.ind, 7, 7),
#                         nbe.genes = nbe.genes)
#  

## ----reduce.signal, eval = FALSE-----------------------------------------
#  # reduced signal by half
#  group.id <- substr(colnames(smp.eff.nc), 7, 7)
#  redhalf.smp.eff.nc <- reduce.signal(smp.eff = smp.eff.nc,
#                              group.id = group.id,
#                              group.id.level = c("E", "V"),
#                              reduce.multiplier = 1/2)
#  
#  
#  # extract differential expressed genes
#  smp.eff.nc.psl <- med.sum.pbset(smp.eff.nc)
#  s.e.limma.fit <- limma.pbset(data = smp.eff.nc.psl,
#                           group.id = group.id,
#                           group.id.level = c("E", "V"))
#  de.ind <- s.e.limma.fit$P.Value < 0.01
#  
#  de.gene.name <- rownames(redhalf.smp.eff.nc)[which(de.ind)][2]
#  nonde.gene.name <- rownames(redhalf.smp.eff.nc)[-which(de.ind)][2]
#  
#  smp.eff.nc.de <- smp.eff.nc[de.gene.name, ]
#  redhalf.smp.eff.nc.de <- redhalf.smp.eff.nc[de.gene.name, ]
#  
#  smp.eff.nc.nonde <- smp.eff.nc[nonde.gene.name, ]
#  redhalf.smp.eff.nc.nonde <- redhalf.smp.eff.nc[nonde.gene.name, ]
#  

## ----reduce.signal.fig, echo = FALSE, fig.align = "hold", fig.width = 4, fig.height = 4----
group.id <- substr(colnames(smp.eff.nc), 7, 7)

# plot
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))
# for DE genes, the difference has been shrunken in half...
boxplot(smp.eff.nc.de ~ group.id == "E", main = "DE")
boxplot(redhalf.smp.eff.nc.de ~ group.id == "E", main = "DE")
# for non-DE genes, the difference has not been changed
boxplot(smp.eff.nc.nonde ~ group.id == "E", main = "Non-DE")
boxplot(redhalf.smp.eff.nc.nonde ~ group.id == "E", main = "Non-DE")


## ----amplify.ary.eff, fig.align = "hold", fig.width = 7, fig.height = 7----
ary.eff.nc.tr <- ary.eff.nc[, c(1:64, 129:192)]

# shift
ary.eff.nc.tr.shift <- amplify.ary.eff(ary.eff = ary.eff.nc.tr,
       amplify.ary.id = colnames(ary.eff.nc.tr)[1:64],
       amplify.level = 2, type = "shift")

# scale 1
ary.eff.nc.tr.scale1 <- amplify.ary.eff(ary.eff = ary.eff.nc.tr,
       amplify.ary.id = colnames(ary.eff.nc.tr)[1:64],
       amplify.level = 2, type = "scale1")

# scale 2
amplify.ary.id <- list(1:40, 41:64, (129:160) - 64, (161:192) - 64)
for(i in 1:length(amplify.ary.id)) 
  amplify.ary.id[[i]] <- colnames(ary.eff.nc.tr)[amplify.ary.id[[i]]]
amplify.level <- c(1.2, 1.3, 1/3, 2/3)

ary.eff.nc.tr.scale2 <- amplify.ary.eff(ary.eff = ary.eff.nc.tr,
       amplify.ary.id = amplify.ary.id,
       amplify.level = amplify.level,
       type = "scale2")

par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))
rng <- range(ary.eff.nc.tr, ary.eff.nc.tr.shift, 
             ary.eff.nc.tr.scale1, ary.eff.nc.tr.scale2)
boxplot(ary.eff.nc.tr, main = "Original",
        ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
boxplot(ary.eff.nc.tr.shift, main = "Shift",
        ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
boxplot(ary.eff.nc.tr.scale1, main = "Scaling 1",
        ylim = rng, pch = 20, cex = 0.2, xaxt = "n")
boxplot(ary.eff.nc.tr.scale2, main = "Scaling 2",
        ylim = rng, pch = 20, cex = 0.2, xaxt = "n")


