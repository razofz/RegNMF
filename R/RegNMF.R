#' @param
#' E: Expression
#' O: Openness
#' Symbol: Symbol Name
#'
"RegNMF" <- function(E,
                     O,
                     symbol,
                     peak_name,
                     symbol_location,
                     peak_location,
                     feature_cut_perc = 0.01,
                     corr_cut_k = 100000,
                     cell_num_smooth = sqrt(dim(E)[2]),
                     core = 8) {
  UseMethod("RegNMF")
}
