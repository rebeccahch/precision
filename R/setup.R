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

'
