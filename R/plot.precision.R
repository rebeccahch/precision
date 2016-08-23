# df.int <- data.frame(precision.results$error_store[1, ])
# df.ext.uh <- data.frame(precision.results$error_store[2, ])
# df.ext.sim.nuh <- data.frame(precision.results$error_store[3, ])
#
# ylim <- c(0, 0.4)
# mytitle <- ""
# design.tr.list = c("CC+", "CC+", "STR", "STR");
# design.te.list = c("CC+", "STR", "CC+", "STR");
# norm.list = c("QN", "NN");
# class.list = c("PAM");
# valid.list = c("int", "ext.uh", "ext.sim.nuh");
# design.tr.list2 <- tolower(gsub("[-]", "2", gsub("[+]", "1", design.tr.list)))
# design.te.list2 <- tolower(gsub("[-]", "2", gsub("[+]", "1", design.te.list)))
#
#
# comb.design <- paste(design.tr.list2, design.te.list2, sep = ".")
# (design.order <-  paste(rep(rep(comb.design, each = length(norm.list)), each = length(class.list)),
#                         rep(class.list, each = length(design.tr.list2)*length(norm.list)),
#                         rep(norm.list, length(design.tr.list2)*length(class.list)),
#                         sep = "."))
#
# boxplot(df.ext.uh[, design.order],
#         las = 2, cex.axis = 0.8,
#         ylab = "Misclassification error rate",
#         xlim = c(0.4, ncol(df.ext.uh[, design.order]) + 0.2),
#         ylim = ylim, boxwex = 0.35,
#         col = adjustcolor("skyblue", alpha.f = 0.6),
#         border = "deepskyblue3", xaxt = "n", pch = 4,
#         main = mytitle)
#
# boxplot(df.int[, design.order],
#         at = 1:ncol(df.int[, design.order]),
#         las = 2,
#         add = TRUE, border = "red", col = NA, boxwex = 0.2,
#         xaxt = "n", yaxt = "n", pch = 2)
# grid(nx = NA, ny = NULL, lwd = 2)
#
# labs <- gsub("STR", "Stratification",
#              gsub("BLK", "Blocking",
#                   gsub("PC2", "Partial Confounding 2",
#                        gsub("PC1", "Partial Confounding 1",
#                             gsub("CC2", "Complete Confounding 2",
#                                  gsub("CC1", "Complete Confounding 1",
#                                            gsub("LASSO.", "",
#                                                 gsub("PAM.", "",
#                                                      gsub("NN", ", No Normalization",
#                                                           gsub("QN", ", Quantile Normalization",
#                                                                toupper(design.order)))))))))))
# labs <- gsub("[.]", " -\n", gsub(".,", ",\n", labs))
#
# axis(1, cex.axis = 0.6, at = 1:length(labs),
#      col = NA,
#      labels = labs)
#
# legend("topleft",
#        legend = c("Cross Validation",
#                   "External Validation"),
#        border = c("red", "deepskyblue3"),
#        bg = "white",
#        fill = c(NA, adjustcolor("skyblue", alpha.f = 0.6)))
#
