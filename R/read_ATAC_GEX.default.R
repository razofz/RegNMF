read_ATAC_GEX.default <- function(foldername) {
  if (typeof(foldername) != "character") {
    stop("foldername should be character\n")
  } else {
    filename <- paste0(foldername, "matrix.mtx")
    a <- read.table(filename, sep = "", header = T, comment.char = "%")
    filename <- paste0(foldername, "features.tsv")
    C <- read.table(filename, sep = "\t")
    chr <- unique(C[, 4])
    chr <- chr[grep("^chr", chr)]
    is_atac <- match(C[, 3], "Peaks")
    is_atac[is.na(is_atac)] <- 0
    is_atac <- as.logical(is_atac)
    peak_name <- C[is_atac, 2]
    symbol <- C[!is_atac, 2]
    features <- C[, 2]

    chr_idx <- match(C[, 4], chr)
    chr_idx[is.na(chr_idx)] <- 0

    symbol_location <- matrix(chr_idx[!is_atac],
      nrow = length(chr_idx[!is_atac])
    )
    symbol_location <- cbind(symbol_location, C[!is_atac, 5])
    peak_name_location <- matrix(chr_idx[is_atac],
      nrow = length(chr_idx[is_atac])
    )
    peak_name_location <- cbind(peak_name_location, C[is_atac, 5])
    filename <- paste0(foldername, "barcodes.tsv")
    barcode <- read.table(filename, sep = "\t")

    ## rna
    f <- match(features, symbol)
    d <- f
    d[is.na(d)] <- 0
    d <- as.logical(d)

    tmp <- match(a[, 1], which(d == TRUE))
    tmp[is.na(tmp)] <- 0
    tmp <- as.logical(tmp)
    a_rna <- a[tmp, ]
    a_rna[, 1] <- f[a_rna[, 1]]
    E <- sparseMatrix(a_rna[, 1], a_rna[, 2],
      x = log2(1 + a_rna[, 3]),
      dims = c(length(symbol), length(barcode[, 1]))
    )

    ## atac
    f <- match(features, peak_name)
    d <- f
    d[is.na(d)] <- 0
    d <- as.logical(d)

    tmp <- match(a[, 1], which(d == TRUE))
    tmp[is.na(tmp)] <- 0
    tmp <- as.logical(tmp)
    a_atac <- a[tmp, ]
    a_atac[, 1] <- f[a_atac[, 1]]
    O <- sparseMatrix(a_atac[, 1], a_atac[, 2],
      x = log10(1 + a_atac[, 3]),
      dims = c(length(peak_name), length(barcode[, 1]))
    )

    ##
    idx <- symbol_location[, 1] > 0
    symbol <- symbol[idx]
    symbol_location <- symbol_location[idx, ]
    E <- E[idx, ]
    idx <- peak_name_location[, 1] > 0
    peak_name <- peak_name[idx]
    peak_name_location <- peak_name_location[idx, ]
    O <- O[idx, ]
    ## non-zero 10 cell
    idx <- rowSums(as.matrix(E) > 0) > 0
    symbol <- symbol[idx]
    symbol_location <- symbol_location[idx, ]
    E <- E[idx, ]
    idx <- rowSums(as.matrix(O) > 0) > 10
    peak_name <- peak_name[idx]
    peak_name_location <- peak_name_location[idx, ]
    O <- O[idx, ]

    return(list(
      E = E,
      O = O,
      PeakName = peak_name,
      barcode = barcode,
      Symbol_location = symbol_location,
      Symbol = symbol,
      Peak_location = peak_name_location,
      chr = chr
    ))
  }
}
