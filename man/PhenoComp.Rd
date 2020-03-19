\name{PhenoComp}
\alias{PhenoComp}
\title{Identification of population-level differentially expressed genes in one-phenotype data}
\usage{
PhenoComp(expdata,label,gene,freq,method,freq1,outfile1,outfile2)
}
\arguments{
  \item{expdata}{a (non-empty) numeric gene expression matrix with both disease and control samples.}
  \item{label}{a (non-empty) numeric vector of data values where ’0’ represents control sample label and ’1’ reptesents disease sample(default).The length of label must be the same as the number of columns in the expdata.}
  \item{gene}{a (non-empty) numeric vector of Entrez gene IDs. The length of gene must be the same as the number of rows in the expdata}
  \item{freq}{the criteria for identifying stable gene pairs in control samples. The default setting of freq is 0.99.}
  \item{method}{Method determines how to estimate the p_up and p_down. Method=1: the p_up and p_down were estimated as the median values of the frequencies of up-regulated and down-regulated genes for individual disease samples.Method=2: the p_up and p_down were estimated as the mean values of the frequencies of up-regulated and down-regulated genes for individual disease samples.
}
  \item{freq1}{the threshold of FDR for identifying population-level differentially expressed genes.}
  \item{outfile1}{The file name used to save the identified population-level up-regulation genes.}
  \item{outfile2}{The file name used to save the identified population-level down-regulation genes.}
}
\description{
PhenoComp is an algorithm to identify population-level differential genes in one-phenotype data. This algorithm is based on RankComp, an algorithm used to identify individual-level differentially expressed genes in each sample.
}
\examples{
PhenoComp(expdata,label,gene,0.99,1,0.05,"gene_up.txt","gene_down.txt")
}
