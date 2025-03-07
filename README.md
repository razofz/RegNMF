# scREG package for singel cell multiome gene expression and chromatin accessibility data
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5676303.svg)](https://zenodo.org/record/5676303#.YY2FWGDMJaQ)

## install package

```R
library(withr)
setRepositories(ind=1:3)
devtools::install_github("Durenlab/RegNMF",ref="main")
```

## Requirements

### system

* macs2(See the detail from [github](https://github.com/macs3-project/MACS))

Check these path that will be used in our program.

```bash
#Using "which" to check the path
which macs2
which awk
```

### dataset

1. Go to [database](https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets) choose one of datasets.
2. You may have to input your infomation for downloading the dataset.
3. Download "Filtered feature barcode matrix MEX (DIR)" and "ATAC Per fragment information file (TSV.GZ)".
4. Unzip them

```bash
#On linux you can use these commands to unzip
    #Unzip Filtered feature barcode matrix MEX
    tar -zxvf [XXX]filtered_feature_bc_matrix.tar.gz
    #There will be a folder named "filtered_feature_bc_matrix/" contain "barcodes.tsv.gz", "matrix.mtx.gz", "features.tsv.gz". Unzip them
    gunzip filtered_feature_bc_matrix/*.gz

    #Unzip ATAC Per fragment information file
    gunzip [XXX]atac_fragments.tsv.gz
```

We will use these data in our program.

## argument

* in_foldername (Character) : Path of unziped Filtered feature barcode matrix MEX("filtered_feature_bc_matrix")
* out_foldername (Character) : Path of folder contain result.
* fragment (Character) : Path of unziped ATAC Per fragment information file("XXXatac_fragments.tsv")
* macs2path (Character) : Path of macs2
* awkpath (Character) : Path of bedtools
* chr (Character) : Which chromatin you want to see in the result(ex. "chr16").
* from (int), to (int) : Which region of the chromasome in the result.
* core (int) : How many core you want to use. You can use `detectCores()` function to check how many core you can use in R.
* width (int), height (int) : Figure size of result.

## simple usage

You can use demo function for clustering.  
See the function detail at ./man/demo.Rd

```R
demo(in_foldername,out_foldername,fragment,macs2path,awktoolspath,core)
```

Ensure there are files named "matrix.mtx", "features.tsv", "barcodes.tsv" in the input folder.

## Using the individual functions  

`callpeak()` in the fourth step requiers MACS2 for peak calling, so you may not able to run the fourth step if you have not installed MACS2.

### First step

Use "read_ATAC_GEX" for loading data.

```R
element=read_ATAC_GEX(in_foldername)
```

### Second step

Use "RegNMF" for cross-modalities dimension reduction; use "clustering" for subpopulation identification.

```R
W123H=RegNMF(E=element$E, 
             O=element$O, 
             Symbol=element$Symbol, 
             PeakName=element$PeakName, 
             Symbol_location=element$Symbol_location, 
             Peak_location=element$Peak_location)

ans=clustering(W123H$H)
```

ans$plot is a figure of tsne.

### Third step

Use "SplitGroup" to predict raw subpoplation-specific cis-regulatory networks. These networks could be further refined in the next step.

```R
groupName=SplitGroup(foldername=out_foldername,
                     barcord=element$barcode[,1],
                     W3=W123H$W3,
                     H=W123H$H,
                     Reg_symbol_name=W123H$Reg_gene_name,
                     Reg_peak_name=W123H$Reg_peak_name,
                     cluster=ans$S[1,])
```

In this function we'll make a file shows pair  of barcords and clusters and a folder contains predicted regulations in each clusters.

### Fourth step

Use "callpeak" to call peak in each clusters and refine the subpopulation-specific cis-regulatory netyworks.

```R
visual_need=callpeak(outfolder=out_foldername,
                     fragment=fragment,
                     barcord_cluster_whole=groupName["barcordFileName"],
                     oldRegFolder=groupName["RegFolderName"],
                     macs2path=macs2path,
                     awkpath=awkpath,
                     cluster=unique(ans$S[1,]),
                     clusterL=length(ans$S[1,]))
```

In this function we'll make three folders, which named "barcord_cluster", "peak_cluster"and "RE_cluster" , contain barcords, peaks, regulations infomation by each clusters.

### Fifth step

Use "Visualization" to run interactive visualization function that takes the genomics region as input and plots the genes, peaks, and interactions in the given range. It includes the genes, the raw peaks from all cells (before clustering), peaks of each cluster from MACS2, and the predicted (refined) peak-gene association in each cluster . The result will be output as "result.pdf"

```R
Visualization(wholef=in_foldername,
              peakf=visual_need["peak_clusterF"],
              regf=visual_need["RE_clusterF"],
              chr=chr,
              from=from,
              to=to,
              clusterlist=clusterlist,
              width=width,
              height=height)
```
