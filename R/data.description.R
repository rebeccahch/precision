#' @title The uniformly-handled probe-level dataset,
#' 10 probes for each unique probe
#'
#' @description A five percent random subset of the uniformly-handled
#' probe-level dataset, 10 probes per each unique probe.
#' The expressions are on a log2 scale without background adjustmnet.
#' This dataset consists of 181 unique probes, of which
#' 6 are negatively biological control probes from Agilent array platform:
#' "NC2_00079215", "NC1_00000215", "NC1_00000197", "NC2_00122731",
#' "NC2_00092197", and "NC2_00106057".
#' The sample IDs (the column names) ending with "E" or "V" are used to indicate
#' whether a sample is endometrial or ovarian tumor sample. There are
#' 96 endometrial and 96 ovarian tumor samples.
#' @format A data matrix with 1810 rows (probes) and 192 columns (samples).
#' @keywords example.data

"uhdata.pl"

#' @title The nonuniformly-handled probe-level dataset,
#' 10 probes for each unique probe
#'
#' @description A five percent random subset of the nonuniformly-handled
#' probe-level dataset, 10 probes per each unique probe.
#' The expressions are on a log2 scale without background adjustmnet.
#' This dataset consists of 181 unique probes, of which
#' 6 are negatively biological control probes from Agilent array platform:
#' "NC2_00079215", "NC1_00000215", "NC1_00000197", "NC2_00122731",
#' "NC2_00092197", and "NC2_00106057".
#' The sample IDs (the column names) ending with "E" or "V" are used to indicate
#' whether a sample is endometrial or ovarian tumor sample. There are
#' 96 endometrial and 96 ovarian tumor samples.
#' @format A data matrix with 1810 rows (probes) and 192 columns (samples).
#' @keywords example.data

"nuhdata.pl"

