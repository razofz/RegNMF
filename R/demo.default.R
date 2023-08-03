# demo.default <- function(in_foldername,out_foldername,fragment,macs2path,bedtoolspath,chr,from,to,core,width=6,height=4){
demo.default <- function(in_foldername,
                         out_foldername,
                         fragment,
                         macs2path,
                         awkpath,
                         core) {
  element <- read_ATAC_GEX.default(in_foldername)
  W123H <- RegNMF.default(
    E = element$E,
    O = element$O,
    symbol = element$symbol,
    peak_name = element$peak_name,
    symbol_location = element$symbol_location,
    peak_location = element$peak_location
  )

  ans <- clustering.default(W123H$H)
  group_name <- SplitGroup.default(
    foldername = out_foldername,
    barcode = element$barcode[, 1],
    W3 = W123H$W3,
    H = W123H$H,
    reg_symbol_name = W123H$reg_gene_name,
    reg_peak_name = W123H$reg_peak_name,
    cluster = ans$S[1, ]
  )

  visual_need <- callpeak.default(
    outfolder = out_foldername,
    fragment = fragment,
    barcode_cluster_whole = group_name["barcode_filename"],
    old_reg_folder = group_name["reg_folder_name"],
    macs2path = macs2path,
    awkpath = awkpath,
    cluster = unique(ans$S[1, ]),
    clusterL = length(ans$S[1, ])
  )

  #  clusterlist=unique(ans$S[1,])

  #  Visualization.default(wholef=in_foldername,
  #                        peakf=visual_need["peak_clusterF"],
  #                        regf=visual_need["RE_clusterF"],
  #                        chr=chr,
  #                        from=from,
  ##                        to=to,
  #                        clusterlist=clusterlist,
  #                        width=width,height=height)
  return(ans)
}
