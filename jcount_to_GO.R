# Load libraries
library(tidyverse)
library(clusterProfiler)  # For GSEA
library('org.Sc.sgd.db')  # Database for Gene Onthology
library(gprofiler2)  # For Over Representation Analysis
library(multimode) 

# Read table made with https://www.yeastgenome.org/
# Contains systematic and some standard gene names
to <- read_tsv("tight_orfs.tsv") %>%
  transmute(
    sys_name = `Gene > Systematic Name`, 
    std_name = `Gene > Standard Name`
  )

# Read table with other synonyms of gene names. https://www.yeastgenome.org/
yeast_genes <- read_tsv('yeast_gene_names.txt')

# Function to convert all gene names either to systematic or standard with to_std=T option
convert_names <- function(vals, to_std = F){
  if (length(vals) == 0){return(vals)}
  l <- tibble(syn = vals, id = (1:length(vals))) %>%
    separate_longer_delim(syn, '|') %>%
    left_join(yeast_genes, by = 'syn')
  
  if (to_std) {
    return(
      l %>%
        mutate(res = ifelse(is.na(std_name), syn, std_name)) %>%
        group_by(id) %>%
        summarise(res = str_flatten(res %>% unique, collapse = '|')) %>%
        pull(res)
    )
  }
  return(
    l %>%
      mutate(res = ifelse(is.na(sys_name), syn, sys_name)) %>%
      group_by(id) %>%
      summarise(res = str_flatten(res %>% unique, collapse = '|')) %>%
      pull(res)
  )
}

#three following functions are presented as example of "cut_significant" function in ORA (look below)
llk <- function(params, data) {
  mu <- params[1]
  sigma <- params[2]
  # calculate log P
  ll <- sum(dnorm(data, mean = mu, sd = sigma, log = TRUE))
  return(-ll)
}


opt_llh<-function(x){
  m<-0
  optim(par=c(m,1), llk, data=x)$par
  
}

######### function to cut significant results
cut_significant <- function(df, col, log = T){
  vals <- df %>% pull({{col}})
  if (log) {vals <- log(vals)}
  mu_sigma <- opt_llh(vals)
  quantile <- qnorm(0.95, mean = mu_sigma[1], sd = mu_sigma[2])
  df %>% filter(vals > quantile)
}

# Paste0 alias
"%+%" <- function(...){
  paste0(...)
}


jcounts <- read_csv("jcounts_exp.csv")

# Normalise the counts
jc <- jcounts %>%
  filter(if_any(where(is.numeric), ~.x!=0))
n = jc %>% dplyr::select(!name) %>% colSums %>% mean
jc <- jc %>% mutate(across(!name, ~ .x / sum(.x) * n))


jc %>%
  filter(d1 > 0) %>%
  mutate(d = log2(d1)) %>%
  mutate(delta = max(d) - (max(d) - locmodes(d)$locations) * 2) %>%
  filter(d > delta) %>%
  pull(d1) %>%
  min() -> min_count

jc%>%
  filter(d1 >= min_count)%>%
  mutate(across(!name, ~ .x + 1))%>%
  mutate(across(!name, ~ .x / d1))%>%
  mutate(across(!name, ~ log2(.x))) -> to_go

to_go$g1 %>% hist(breaks = 1000)


################
##### GSEA
################

# Initialize GO_norm tibble
GO_norm <- tibble(
  pvalue = numeric(), 
  Description = character(), 
  condition = character(), 
  NES = numeric(), 
  geneID = character()
)

# Run GSEA for various conditions
for (i in c('g1_d1', 'g2_d2', 'd2_d1', 'g2_g1')) {
  
  # Define condition
  condition <- i
  print(paste("GSEA", condition))
  
  # Split the condition into t and u
  t = str_split_1(condition, '_')[1]
  u = str_split_1(condition, '_')[2]
  
  # Perform mutations and select columns
  to_go %>%
    dplyr::select(name, starts_with(t) | starts_with(u)) %>%
    mutate(log2FoldChange = !!sym(t) - !!sym(u)) %>% 
    mutate(gene = name) %>% 
    dplyr::select(gene, log2FoldChange) -> res
  
  # Consider cutting conditions
  cut <- res %>% filter(!str_detect(gene, '\\|'))
  
  # Process gene list
  gene_list <- cut$log2FoldChange
  names(gene_list) <- cut$gene %>% convert_names
  gene_list <- sort(gene_list, decreasing = T)
  
  # Run GSEA
  norm_gsea <- clusterProfiler::gseGO(
    geneList = gene_list,
    OrgDb = org.Sc.sgd.db,
    keyType = "ENSEMBL",
    ont = "BP"
  )
  
  # Add result to GO_norm if "core_enrichment" is present
  if ( "core_enrichment" %in% names(norm_gsea@result)) {
    gsea1 <- norm_gsea@result%>%
      as_tibble()%>%
      group_by(core_enrichment)%>%
      filter(setSize == min(setSize))%>%
      ungroup()
    
    gsea1$core_enrichment%>%str_split('\\/') -> genesets
    m <- sapply(genesets, 
                function(j){
                  sapply(genesets, 
                         function(i){
                           res <- sum(i %in% j) == length(i)
                           if((res == T) & (res == (sum(j %in% i) == length(j)))){
                             res <- FALSE
                           }
                           res })
                })#true means j contains i
    colSums(m) -> par
    res <- norm_gsea%>%
      filter(ID %in% filter(gsea1, par == 0)$ID)    
    
    gene_clust <- res@result$core_enrichment %>%
      map(~convert_names(str_split(., '\\/') %>% unlist(), to_std = T)) %>%
      map(~str_flatten(., collapse = '/')) %>% as.character()
    
    GO_norm <- res@result %>%
      dplyr::select(Description, NES, pvalue) %>%
      mutate(condition = condition, geneID = gene_clust) %>%
      add_row(GO_norm,.)
  }
}

# Write the results to a TSV file
GO_norm %>% 
  arrange(condition, NES, pvalue) %>%
  write_tsv('./GSEA_exp.txt')

# Read the TSV file
GSEA_exp_to_split <- read_tsv('./GSEA_exp.txt')

# Get unique conditions
files <- GSEA_exp_to_split$condition %>% unique

# Write each condition to a separate file
map(files, ~ write_tsv(GSEA_exp_to_split %>% 
                         filter(str_detect(condition, .x)),
                       './GSE_exp_' %+% .x %+% '.txt'))

################
##### ORA
################
GO_report <- tibble(p_value = numeric(), 
                    term_name = character(),  
                    condition = character(), 
                    intersection = character(), 
                    overrep = logical())

for (i in c('g1_d1', 'g2_d2', 'd2_d1', 'g2_g1')){
  condition <- i
  print(paste("ORA",condition))
  t = str_split_1(condition, '_')[1]
  u = str_split_1(condition, '_')[2]
  to_go %>%
    dplyr::select(name, starts_with(t) | starts_with(u)) %>%
    mutate(log2FoldChange = !!sym(t) - !!sym(u)) %>% 
    mutate(gene = name) %>% 
    dplyr::select(gene, log2FoldChange) -> res

res %>% cut_significant(log2FoldChange, F) -> cut

geneList <- cut$gene%>%convert_names()
bgList <- res$gene%>%convert_names()

GOres <- gost(geneList, 
              organism = 'scerevisiae', 
              ordered_query = T, 
              measure_underrepresentation = F, 
              user_threshold = 0.05, 
              significant = F, 
              correction_method = 'fdr', 
              custom_bg = bgList, 
              sources = c('GO:BP', 'KEGG'), 
              domain_scope = 'custom_annotated', 
              evcodes = TRUE)
GO_report <- GOres$res%>%
  dplyr::select(p_value, term_name, intersection)%>%
  as_tibble()%>%
  mutate(condition = condition, 
         intersection = str_replace_all(intersection, ',', '|')%>%
           convert_names(to_std = T))%>%
  arrange(p_value)%>%
  mutate(overrep = T)%>%add_row(GO_report,.)


res%>%
  mutate(log2FoldChange = - log2FoldChange)%>%
  cut_significant(log2FoldChange, F)->cut

geneList <- cut$gene %>% convert_names()

GOres <- gost(geneList, 
              organism = 'scerevisiae', 
              ordered_query = T, 
              measure_underrepresentation = F, 
              user_threshold = 0.05, 
              significant = F, 
              correction_method = 'fdr', 
              custom_bg = bgList, 
              sources = c('GO:BP', 'KEGG'), 
              domain_scope = 'custom_annotated', 
              evcodes = TRUE)
GO_report <- GOres$res%>%
  dplyr::select(p_value, term_name, intersection)%>%
  as_tibble()%>%
  mutate(condition = condition,
         intersection = str_replace_all(intersection, ',', '|')%>%
           convert_names(to_std = T))%>%
  arrange(p_value)%>%
  mutate(overrep = F)%>%add_row(GO_report,.)
}

write_tsv(GO_report%>%arrange(p_value), './exp_GO_report.txt')
