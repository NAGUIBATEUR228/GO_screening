# Gene onthology of yeast deletion collection BarSeq analysis

Yeast deletion (knockout) collection - **YKOC** - is used for genetic screenings, because each strain has a deletion with a unique 20-nucleotide barcode flanked with adapters for primers. Barcode suquencing (BarSeq) allows to estimate relative ammount of each deletion strain in different conditions. But this information is not full enough for choosing any strain for following assay.
Here we introduce a pipeline for Gene Onthology analysis of BarSeq experiment data.

The code is written in R 4.3.1 language


```r
library(tidyverse)
library("ggrepel")
library('patchwork')
library(clusterProfiler)#for GSEA
library('org.Sc.sgd.db')
library(gprofiler2)#for Fischer enrichment
```
