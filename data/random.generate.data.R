set.seed(1001)

uh <- matrix(rnorm(n = length(nuhdata.pl), mean = 0, sd = 0.3),
              ncol = ncol(nuhdata.pl), nrow = nrow(nuhdata.pl))

nuh <- uh + matrix(rnorm(n = length(nuhdata.pl), mean = 0, sd = 1),
                   ncol = ncol(nuhdata.pl), nrow = nrow(nuhdata.pl))

cn <- paste0(substr(colnames(nuhdata.pl), 1, 2),
             sample(seq(4000, 6000), 192, replace = FALSE),
             substr(colnames(nuhdata.pl), 7, 7))

colnames(uh) <- sample(cn)
colnames(nuh) <- cn
rownames(uh) <- rownames(nuhdata.pl)
rownames(nuh) <- rownames(nuhdata.pl)

boxplot(uh, las = 2, pch = 20, cex = 0.5, cex.lab = 0.3)
boxplot(nuh, las = 2, pch = 20, cex = 0.5, cex.lab = 0.3)

uhdata.pl <- uh
nuhdata.pl <- nuh
