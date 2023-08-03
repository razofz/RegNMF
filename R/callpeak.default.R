callpeak.default <- function(outfolder,
                             fragment,
                             barcode_cluster_whole,
                             old_reg_folder,
                             macs2path,
                             awkpath,
                             cluster,
                             clusterL) {
  if (!str_ends(outfolder, pattern = "/")) {
    outfolder <- str_c(outfolder, "/")
  }
  barcode_clusterF <- paste0(outfolder, "barcode_cluster/")
  if (!dir.exists(barcode_clusterF)) {
    dir.create(barcode_clusterF, recursive = TRUE)
  } else if (length(dir(barcode_clusterF))) {
    file.remove(paste0(barcode_clusterF, dir(barcode_clusterF)))
  }

  RE_clusterF <- paste0(outfolder, "RE_cluster/")

  if (!dir.exists(RE_clusterF)) {
    dir.create(RE_clusterF, recursive = TRUE)
  } else if (length(dir(RE_clusterF))) {
    file.remove(paste0(RE_clusterF, dir(RE_clusterF)))
  }

  peak_clusterF <- paste0(outfolder, "peak_cluster/")
  if (!dir.exists(peak_clusterF)) {
    dir.create(peak_clusterF, recursive = TRUE)
  } else if (length(dir(peak_clusterF))) {
    file.remove(paste0(peak_clusterF, dir(peak_clusterF)))
  }

  # print("hey")
  cmd <- paste0(
    awkpath,
    " -v l=",
    clusterL,
    " -v folder=",
    barcode_clusterF,
    " \'{if(NR<=l){col[$1]=$2} else{if(length(col[$4])!=0)",
    "{print >>folder col[$4]\".bed\"}}}\' ",
    barcode_cluster_whole,
    " ",
    fragment
  )
  system(command = cmd)

  for (i in cluster) {
    cmd <- paste0(
      macs2path,
      " callpeak -t ",
      barcode_clusterF,
      i,
      ".bed -g hs -f BED --nomodel --shift -100 --extsize 200 -n ",
      peak_clusterF,
      i
    )
    system(command = cmd)
    genome <- Seqinfo(genome = NA_character_)
    peakfile <- paste0(peak_clusterF, i, "_peaks.narrowPeak")
    oldREfile <- paste0(old_reg_folder, "reg_cluster", i, ".bed")
    REfile <- paste0(RE_clusterF, i, ".bed")
    gr_a <- import(peakfile, genome = genome)
    gr_b <- import(oldREfile, genome = genome)
    pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE)
    df <- data.frame(
      chr = pairs@first@seqnames,
      start = pairs@first@ranges@start,
      end = (pairs@first@ranges@start + pairs@first@ranges@width),
      TG = pairs@second@elementMetadata@listData$name,
      score = pairs@second@elementMetadata@listData$score
    )
    write.table(df,
      REfile,
      sep = "\t",
      col.names = F,
      row.names = F,
      quote = FALSE
    )
  }

  return(c(RE_clusterF = RE_clusterF, peak_clusterF = peak_clusterF))
}
