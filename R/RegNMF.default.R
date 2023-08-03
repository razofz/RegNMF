RegNMF.default <- function(E,
                           O,
                           symbol,
                           peak_name,
                           symbol_location,
                           peak_location,
                           feature_cut_perc = 0.01,
                           corr_cut_k = 100000,
                           cell_num_smooth = sqrt(dim(E)[2]),
                           core = 8) {
  # K parameter
  K <- 100
  ##### Normalize and filter
  E1 <- normalize_scaling.default(E)
  O1 <- normalize_scaling.default(O)
  numCut <- feature_cut_perc * dim(E1)[2]
  a <- Matrix::rowSums(E1 > 0)
  a1 <- Matrix::rowSums(O1 > 0)
  E11 <- E1[a > numCut, ]
  rm(E1)
  gc()

  O11 <- O1[a1 > numCut, ]
  rm(O1)
  gc()

  symbol_location <- symbol_location[a > numCut, ]
  peak_location <- peak_location[a1 > numCut, ]
  symbol <- symbol[a > numCut]
  peak_name <- peak_name[a1 > numCut]
  rm(a1)
  gc()
  ls()

  ### Making reg matrix
  R <- Fold_RE_TG_multiAdjust.default(E11, O11, symbol_location, peak_location)
  R <- as(R, "sparseMatrix")
  mk <- sort(R[R > 0], decreasing = TRUE)
  corr_cut_k <- min(length(mk), corr_cut_k)
  mk <- mk[1:corr_cut_k]

  corr_cut <- mk[length(mk)]
  a <- Matrix::colSums(R >= corr_cut)
  O11 <- O11[a > 0, ]
  a2 <- Matrix::rowSums(R >= corr_cut)
  E11 <- E11[a2 > 0, ]
  R <- R[a2 > 0, a > 0]
  symbol_location <- symbol_location[a2 > 0, ]
  peak_location <- peak_location[a > 0, ]
  symbol <- symbol[a2 > 0]
  peak_name <- peak_name[a > 0]
  rm(a)
  gc()
  c_mat <- which(((as.matrix(R) >= corr_cut) == 1), arr.ind = TRUE)
  rm(R)
  gc()
  reg <- O11[c_mat[, 2], ] + E11[c_mat[, 1], ]
  reg_gene_location <- symbol_location[c_mat[, 1], ]
  reg_peak_location <- peak_location[c_mat[, 2], ]
  reg_dis <- abs(reg_gene_location[, 2] - reg_peak_location[, 2])
  # distance penalty parameter d_0 defined here
  d0 <- 2 * 10^5
  reg_w <- exp(-reg_dis / d0)
  reg_adj <- reg * reg_w
  O12 <- tfidf_trans.default(O11)
  rm(O11)
  gc()

  reg_info <- cbind(c_mat[, 1], c_mat[, 2])
  reg_info <- cbind(reg_info, reg_dis)

  ### Making W10,W20,W30,H0, which are used in initialize lambda
  W10 <- matrix(runif(dim(O12)[1] * K, min = 0, max = 1), dim(O12)[1], K)
  W20 <- matrix(runif(dim(E11)[1] * K, min = 0, max = 1), dim(E11)[1], K)
  W30 <- matrix(
    runif(dim(reg_adj)[1] * K, min = 0, max = 1),
    dim(reg_adj)[1],
    K
  )
  H0 <- matrix(runif(dim(O12)[2] * K), K, dim(O12)[2])
  wh1 <- CNmf(as.matrix(O12), 100, 100, W10, H0, core)
  wh2 <- CNmf(as.matrix(E11), 100, 100, W20, H0, core)
  wh3 <- CNmf(as.matrix(reg_adj), 100, 100, W30, H0, core)
  wh1$W <- t(t(wh1$W) * sqrt(rowSums(wh1$H * wh1$H)))
  wh2$W <- t(t(wh2$W) * sqrt(rowSums(wh2$H * wh2$H)))
  wh3$W <- t(t(wh3$W) * sqrt(rowSums(wh3$H * wh3$H)))
  lambda <- defaultpar_CoupledNMF.default(
    O12, wh1$W, E11, wh2$W, reg_adj,
    wh3$W, 2, 1
  )

  rm(wh1)
  rm(wh2)
  rm(wh3)
  gc()

  W10 <- matrix(runif(dim(O12)[1] * K, min = 0, max = 1), dim(O12)[1], K)
  W20 <- matrix(runif(dim(E11)[1] * K, min = 0, max = 1), dim(E11)[1], K)
  W30 <- matrix(
    runif(dim(reg_adj)[1] * K, min = 0, max = 1), dim(reg_adj)[1],
    K
  )
  H0 <- matrix(runif(dim(O12)[2] * K), K, dim(O12)[2])
  W123H <- CPPNMF_cluster_joint_cross_domain_try(
    as.matrix(O12), as.matrix(E11),
    as.matrix(reg_adj), K, 300,
    lambda[[1]], lambda[[2]], W10,
    W20, W30, H0, c_mat[, 1], c_mat[, 2],
    reg_w, core
  )

  rm(W10)
  rm(W20)
  rm(W30)
  rm(O12)
  rm(E11)
  gc()

  W123H[["reg_gene_name"]] <- symbol[c_mat[, 1]]
  W123H[["reg_peak_name"]] <- peak_name[c_mat[, 2]]
  W123H[["gene_name"]] <- symbol
  W123H[["peak_name"]] <- peak_name
  W123H[["reg_adj"]] <- reg_adj

  return(W123H)
}
