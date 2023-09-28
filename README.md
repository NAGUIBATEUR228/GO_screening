# Gene onthology of yeast deletion collection BarSeq analysis

Yeast deletion (knockout) collection - **YKOC** - is used for genetic screenings, because each strain has a deletion with a unique 20-nucleotide barcode flanked with adapters for primers. Barcode suquencing (BarSeq) allows to estimate relative ammount of each deletion strain in different conditions. But this information is not full enough for choosing any strain for following assay.  
Here we introduce a pipeline for Gene Onthology analysis of BarSeq experiment data.

### Repository overview

`barcode_analysis.R` - script that transforms fastq files to count tables  
it requires `puddu.txt` and `ref_nm.txt` and produces also  `reference.txt`:  
`puddu.txt` - table based on Puddu et al. whole genome sequencing of YKOC  
`ref_nm.txt` - a table with original reference barcodes  
`reference.txt` - summary reference barcode table based on `puddu.txt` and `ref_nm.txt`  

`barocde_dm.py` - script that extracts barcodes as `barcode_analysis.R` from reads with inaccurate adapter sequence  
`barcode_blast.py` - script that matches mutant barcodes using BLAST command line program  
> to install blastn on local computer download [here](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) for example `ncbi-blast-2.14.1+-win64.exe` and add bin/blastn to the PATH.  

`exp_to_jcounts.R` - script that joines experimental count tables to ORf vs expriment count table  
`jcount_to_GO.R` -  script that performs GO analysis of `jcounts_exp.txt` table  
both these scripts require `tight_orfs.tsv` and `yeast_gene_names.txt` - the tables with standard and systematic gene names.

`GSEA_exp.txt` - a table with all GO terms found in comparisons
|pvalue|Description|condition|overrep|geneID|
|--|--|--|--|--|
|1e-10|translation|d2_d1|FALSE|ZUO1/RPL21A/LSO2...|

***pvalue*** - corrected pvalue of a GO term  
***Description*** - description of a GO term  
***condition*** - describes which columns of jcounts are compared (treated_untreated)  
***overrep*** - TRUE means that the positive LogFC genes were analysed, i.e. overrepresented in treated. FALSE corresponds to negative LogFC.  
***geneID*** - a list of genes in GO term separated by slash

## [jcount_to_GO.R](https://github.com/NAGUIBATEUR228/GO_screening/blob/main/jcount_to_GO.R)

The code is written in `R 4.3.1` language

```r
library(tidyverse)
library(clusterProfiler) #for GSEA
library('org.Sc.sgd.db') #database for Gene Onthology
library(gprofiler2) #for Over Representation Analysis
```

After writing aiding functions and creating jcounts table the analysis itself starts

Reading file with joined count data
```r
setwd("C:/path/exp")
jcounts <- read_csv("jcounts_exp.csv")
```
|name|d1|d2|g1|g2|
|:--:|:-:|:-:|:-:|:-:|
|AAC1|5857|4949|3399|4859|
|AAC3|2182|7250|1580|12191|
|AAD3|4121|746|6555|1143|

Normalisation of the counts equalizing library size
```r
jc <- jcounts%>%
  filter(if_any(where(is.numeric), ~.x!=0))
n = jc%>%dplyr::select(!name) %>% colSums %>% mean
jc <- jc %>% mutate(across(!name, ~ .x / sum(.x) * n))
```

Evaluating LogFoldChange in comparison to original pool
```r
jc%>%
  filter(d1 > 0)%>%
  mutate(across(!name, ~ .x + 1))%>%
  mutate(across(!name, ~ .x / d1))%>%
  mutate(across(!name, ~ log2(.x))) -> to_go
```

The table for results of Gene Set Enrichment Analysis has columns with p-value of GO-term, name of GO-term, condition i.e. which experiments are compared, overrepresented or not the term in a "treated" experiment, gene IDs in term.
```r
#GSEA
GO_norm <- tibble(pvalue = numeric(),
           Description = character(),
           condition = character(),
           overrep = logical(),
           geneID = character())
```

Subtracting LogFC provides proper LogFC for these experiments in condition.
```r
for (i in c('g1_d1', 'g2_d2', 'd2_d1', 'g2_g1')){
  condition <- i
  t = str_split_1(condition, '_')[1]
  u = str_split_1(condition, '_')[2]
  to_go %>%
    dplyr::select(name, starts_with(t) | starts_with(u)) %>%
    mutate(log2FoldChange = !!sym(t) - !!sym(u)) %>% 
    mutate(gene = name) %>% 
    dplyr::select(gene, log2FoldChange) -> res
```
Cutting a gene list allows to measure separately over- and underrepresented genes. Adding minus to LogFC helps to sort undedrrepresented genes for GSEA. 
```r
 cut <- res %>% filter(!str_detect(gene, '\\|'), log2FoldChange < 0)
  
  gene_list <- -cut$log2FoldChange
  names(gene_list) <- cut$gene%>%convert_names
  gene_list <- sort(gene_list, decreasing = T)
  
  norm_gsea <- clusterProfiler::gseGO(
    geneList = gene_list,
    OrgDb = org.Sc.sgd.db,
    keyType = "ENSEMBL",
    ont = "BP")
```

Look for information on using [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) and [gprofiler2](https://biit.cs.ut.ee/gprofiler/page/r) 



