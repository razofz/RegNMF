SplitGroup.default <- function(foldername,
                               barcode,
                               W3,
                               H,
                               reg_symbol_name,
                               reg_peak_name,
                               cluster) {
  if (!dir.exists(foldername)) {
    dir.create(foldername, recursive = TRUE)
  }
  clustern <- length(unique(cluster))
  barcode_cluster <- data.frame(barcode = barcode, cluster = cluster)
  if (!str_ends(foldername, pattern = "/")) {
    foldername <- str_c(foldername, "/")
  }
  bfilename <- paste0(foldername, "barcode_cluster.bed")
  write.table(barcode_cluster, bfilename,
    sep = "\t", col.names = F, row.names = F, quote = FALSE
  )
  chr <- c()
  peaks <- c()
  peake <- c()
  for (i in 1:length(reg_peak_name)) {
    a <- strsplit(reg_peak_name[i], ":")
    b <- strsplit(a[[1]][2], "-")
    chr[i] <- a[[1]][1]
    peaks[i] <- b[[1]][1]
    peake[i] <- b[[1]][2]
  }
  df <- data.frame(
    chr = chr, peaks = peaks, peake = peake, symbolName = reg_symbol_name
  )

  H_norm <- H / sqrt(rowSums(H * H))
  W3_norm <- t(t(W3) * sqrt(rowSums(H * H)))
  H_w <- matrix(nrow = nrow(H_norm), ncol = clustern)

  for (i in 1:clustern) {
    H_w[, i] <- rowMeans(H_norm[, cluster == i])
  }

  W3_cluster <- W3_norm %*% H_w
  reg_folder_name <- paste0(foldername, "old_reg_cluster/")
  if (!dir.exists(reg_folder_name)) {
    dir.create(reg_folder_name, recursive = TRUE)
  } else if (length(dir(reg_folder_name))) {
    file.remove(paste0(reg_folder_name, dir(reg_folder_name)))
  }
  saveRDS(W3_cluster, str_c(foldername, "W3_cluster.rds"))

  for (i in 1:clustern) {
    topk <- order(W3_cluster[, i], decreasing = TRUE)[1:10000]
    outdf <- df[topk, ]
    outdf$reg <- W3_cluster[topk, i]
    filename <- paste0(reg_folder_name, "reg_cluster", i, ".bed")
    write.table(outdf, filename,
      sep = "\t", col.names = F, row.names = F, quote = FALSE
    )
    topk <- order(W3_cluster[, i], decreasing = TRUE)
    outdf <- df[topk, ]
    outdf$reg <- W3_cluster[topk, i]
    filename <- paste0(reg_folder_name, "reg_cluster", i, "_all.bed")
    write.table(outdf, filename,
      sep = "\t", col.names = F, row.names = F, quote = FALSE
    )
  }

  return(a = c(
    barcode_file_name = bfilename,
    reg_folder_name = reg_folder_name
  ))
}
