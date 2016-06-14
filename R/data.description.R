#' @title The uniformly-handled (benchmark) probe-level dataset, 10 probes for each unique probe
#'
#' @description The uniforly-handled probe-level dataset, with non-control-probe removed,
#' 10 probes per each unique probe, no background adjustment.
#' The expressions are on a log-2 scale.
#'
#' @format A data matrix with 1810 rows (probes) and 192 columns (samples).
#' Its column names ending with "E" or "V" indicate whether a sample belongs to endometrial or ovarian tumor sample.
#' @keywords example.data

"r.data.pl"

#' @title The non-uniformly-handled (test) probe-level dataset, 10 probes for each unique probe
#'
#' @description The non-uniformly-handled probe-level dataset, with non-control-probe removed,
#' 10 probes per each unique probe, no background adjustment. The expressions are on a log-2 scale.
#'
#' @format A data matrix with 1810 rows (probes) and 192 columns (samples).
#' Its column names ending with "E" or "V" indicate whether a sample belongs to endometrial or ovarian tumor sample.
#' @keywords example.data

"non.r.data.pl"

