# library(stats)
# set.seed(1001)
#
# # four random digits for the sample IDs, used to further de-identify the data
# rand.four.digits <- sample(seq(4000, 6000), 192, replace = FALSE)
#
# # add random errors to uhdata.pl
# uhdata.pl.rand <- uhdata.pl.ori +
#   matrix(rnorm(n = length(uhdata.pl.ori), mean = 0, sd = 0.3),
#          ncol = ncol(uhdata.pl.ori), nrow = nrow(uhdata.pl.ori))
#
# cn <- paste0(substr(colnames(uhdata.pl.ori), 1, 2),
#              rand.four.digits,
#              substr(colnames(uhdata.pl.ori), 7, 7))
#
# colnames(uhdata.pl.rand) <- cn
#
# par(mfrow = c(1, 2))
# boxplot(uhdata.pl.ori, las = 2, pch = 20, cex = 0.5, cex.lab = 0.3)
# boxplot(uhdata.pl.rand, las = 2, pch = 20, cex = 0.5, cex.lab = 0.3)
#
#
# # add random errors to nuhdata.pl
# nuhdata.pl.rand <- nuhdata.pl.ori +
#   matrix(rnorm(n = length(nuhdata.pl.ori), mean = 0, sd = 0.3),
#          ncol = ncol(nuhdata.pl.ori), nrow = nrow(nuhdata.pl.ori))
#
# colnames(nuhdata.pl.rand) <- cn[match(colnames(nuhdata.pl.ori),
#                                       colnames(uhdata.pl.ori))]
#
# par(mfrow = c(1, 2))
# boxplot(nuhdata.pl.ori, las = 2, pch = 20, cex = 0.5, cex.lab = 0.3)
# boxplot(nuhdata.pl.rand, las = 2, pch = 20, cex = 0.5, cex.lab = 0.3)
#
#
# uhdata.pl <- uhdata.pl.rand
# nuhdata.pl <- nuhdata.pl.rand
