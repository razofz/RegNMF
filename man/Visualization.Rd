\name{Visualization}
\alias{Visualization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Reg Nonnegative Matrix Factorization
}
\description{
Using the results of callpeak, plot a Figure containing regulation and peak
information.
}
\usage{
  RgeNMF(wholef,
        peakf,
        regf,
        chr,
        from,
        to,
        clusterlist,
        width=6,
        height=4)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{wholef=in_foldername}{Filtered feature barcode matrix file}
\item{peakf}{Folder contain peak information for each cluster}
\item{regf}{Folder contain regulation information for each cluster}
\item{chr}{Which chromosome you want to see in the result(ex. "chr16").}
\item{from,to}{Which region of the chromosome you want to see in the result.}
\item{clusterlist}{Vector of cluster names}
\item{width,height}{Figure size of result.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\author{
Cham, fchang@clemson.edu
}
