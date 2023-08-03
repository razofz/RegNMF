callpeak <- function(outfolder,
                     fragment,
                     barcode_cluster_whole,
                     old_reg_folder,
                     macs2path,
                     awkpath,
                     cluster,
                     clusterL) {
  UseMethod("callpeak")
}
