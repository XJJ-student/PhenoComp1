# PhenoComp
Identification of population-level differentially expressed genes in one-phenotype data

# Install
To install the PhenoComp, install from github using devtools
```
library(devtools)
install_github("XJJ-student/PhenoComp")
```
Or you can download the .ZIP file and and unzip it.
```
install.packages("PhenoComp",repos = NULL,type="source")
#The "PhenoComp" should be combined with the absolute path.
```
# Usage
```
library(PhenoComp)
data(example)
PhenoComp(expdata,label,gene,0.99,1,0.05,"gene_up.txt","gene_down.txt")
```
The example is the gene expression profile of GSE26887 from database Gene Expression Omnibus (GEO)
# Data input
PhenoComp(expdata, label, gene, freq, method, freq1, outfile1, outfile2)
arguments|description
:--|:---
expdata|a (non-empty) numeric gene expression matrix with both disease and control samples.
label|a (non-empty) numeric vector of data values where ’0’ represents control sample label and ’1’ reptesents disease sample(default).The length of label must be the same as the number of columns in the expdata.
gene|a (non-empty) numeric vector of Entrez gene IDs. The length of gene must be the same as the number of rows in the expdata
freq|the criteria for identifying stable gene pairs in control samples. The default setting of freq is 0.99.
method|Method determines how to estimate the p_up and p_down. Method=1: the p_up and p_down were estimated as the median values of the frequencies of up-regulated and down-regulated genes for individual disease samples.Method=2: the p_up and p_down were estimated as the mean values of the frequencies of up-regulated and down-regulated genes for individual disease samples.
freq1|the threshold of FDR for identifying population-level differentially expressed genes.
outfile1|The file name used to save the identified population-level up-regulation genes.
outfile2|The file name used to save the identified population-level down-regulation genes.

# Contact email
Please don't hesitate to address comments/questions/suggestions regarding this R package to:
Jiajing Xie <xiejiajing_fjmu@163.com>; Haidan Yan <Joyan168@126.com>
