Fold_RE_TG_multiAdjust.default <- function(E2,
                                           O2,
                                           symbol_location,
                                           peak_location) {
  E2 <- as.matrix(E2)
  O2 <- as.matrix(O2)
  print(typeof(E2))
  print(typeof(O2))
  print(typeof(symbol_location))
  print(typeof(peak_location))
  P_1 <- Fold_RE_TG_MultiAdjustCore(E2, O2, symbol_location, peak_location)

  return(P_1)
}
