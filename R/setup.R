'
install.packages("devtools")
library(devtools)

install_github("hadley/devtools")
library(devtools)

install.packages("roxygen2")
library(roxygen2)

install_deps()

path <- find.package("PRECISION")
system(paste(shQuote(file.path(R.home("bin"), "R")),"CMD", "Rd2pdf", shQuote(path)))


#** Tutorial **#
library(PRECISION)

## randomized dataset -----
load("/Volumes/BstShared/Biostatistics/Qin/Li-Xuan_Rebecca/Agilent_MiRNA_Array_Data/Workspace/workspace.benchmark.192arrays.final/benchmark.noba.1raw.fin.Rdata")
dim(noba.bench.log2.pl)

# trim to 10 probes per probe-set
noba.bench.log2.pl.p10 <- per.unipbset.truncate(noba.bench.log2.pl)

# summarize probes into each unique probe-set by the median
noba.bench.log2.psl.p10 <- med.sum.pbset(data = noba.bench.log2.pl.p10, num.per.unipbset = 10)

# non-control genes only
bench.psl <- noba.bench.log2.psl.p10[rownames(noba.bench.log2.psl)[probe.type.psl == 0], ]

# run differential expression analysis
group.id <- substr(colnames(bench.psl), 7, 7)
group.id.level <- levels(as.factor(group.id))
limma.fit.bench <- limma.pbset(data = bench.psl, group.id = group.id,
group.id.level = group.id.level)
table(limma.fit.bench$P.Value < 0.01)

# get differentially expressed genes
de.gene <- rownames(limma.fit.bench)[which(limma.fit.bench$P.Value < 0.01)]
non.de.gene <- rownames(limma.fit.bench)[!rownames(limma.fit.bench) %in% de.gene]

# random subset of the data for example dataset
set.seed(1001)
subset.gene <- c(sample(de.gene, size = length(de.gene)*0.05),
sample(non.de.gene, size = length(non.de.gene)*0.05))

# combine negative-biological-control probe data and the random subset data together
ctrl.gene <- rownames(noba.bench.log2.psl)[probe.type.psl == -1][-1]

select.gene <- c(subset.gene, ctrl.gene)
r.data.pl <- noba.bench.log2.pl.p10[rownames(noba.bench.log2.pl.p10) %in% select.gene, ]

## non-randomized dataset -----
load("/Volumes/BstShared/Biostatistics/Qin/Li-Xuan_Rebecca/Agilent_MiRNA_Array_Data/Workspace/workspace.test.190+2arrays.final/test192.noba.1raw.fin.Rdata")
dim(noba.test.log2.pl)

# trim to 10 probes per probe-set
noba.test.log2.pl.p10 <- per.unipbset.truncate(noba.test.log2.pl)
non.r.data.pl <- noba.test.log2.pl.p10[rownames(noba.test.log2.pl.p10) %in% select.gene, ]

dim(non.r.data.pl)



## prestored objects ----
temp <- limma.fit.r.data$P.Value
limma.fit.r.data <- list()
limma.fit.r.data$P.Value <- temp
rm(temp)

temp <- pam.int$mc
temp2 <- pam.pred$mc
pam.int <- pam.pred <- list()
pam.int$mc <- temp
pam.pred$mc <- temp2
rm(temp, temp2)

temp <- lasso.int$mc
temp2 <- lasso.pred$mc
lasso.int <- lasso.pred <- list()
lasso.int$mc <- temp
lasso.pred$mc <- temp2
rm(temp, temp2)

temp <- uni.handled.results$error_store
uni.handled.results <- list()
uni.handled.results$error_store <- temp
rm(temp)

temp <- uni.handled.results.rn$error_store
uni.handled.results.rn <- list()
uni.handled.results.rn$error_store <- temp
rm(temp)

temp <- uni.handled.results.ridge$error_store
uni.handled.results.ridge <- list()
uni.handled.results.ridge$error_store <- temp
rm(temp)

temp <- precision.results$error_store
precision.results <- list()
precision.results$error_store <- temp
rm(temp)

temp <- precision.ruv4.results$error_store
precision.ruv4.results <- list()
precision.ruv4.results$error_store <- temp
rm(temp)

save(limma.fit.r.data, de.ind,
smp.eff.nc.de, redhalf.smp.eff.nc.de,
smp.eff.nc.nonde, redhalf.smp.eff.nc.nonde,
pam.int, pam.pred, lasso.int, lasso.pred,
uni.handled.results, uni.handled.results.rn, uni.handled.results.ridge,
precision.results, precision.ruv4.results,
file = "vignettes/prestore.obj.Rdata")

'
