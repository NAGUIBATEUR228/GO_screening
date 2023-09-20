library(tidyverse)
library("ggrepel")
library('patchwork')
library(clusterProfiler)#for GSEA
library('org.Sc.sgd.db')
library(gprofiler2)#for Fischer enrichment

folder_with_fastq<-'exp'
to<-read_tsv("../tables/tight_orfs.tsv")%>%
  transmute(sys_name=`Gene > Systematic Name`,std_name= `Gene > Standard Name`)

to_join<-function(x){
  y<-x%>%
    transmute(sys_name=Confirmed_deletion,count)%>%
    left_join(to,by=c("sys_name"))%>%
    group_by(sys_name,std_name)%>%
    summarise(count=sum(count))%>%
    ungroup()
  y$count[is.na(y$count)]<-0
  return(y)
}


yeast_genes <- read_tsv('../fin/yeast_gene_names.txt')
convert_names <- function(vals, to_std = F){
  if (length(vals) == 0){return(vals)}
  l <- tibble(syn = vals, id = (1:length(vals)))%>%
    separate_longer_delim(syn, '|')%>%
    left_join(yeast_genes, by = 'syn')
  if (to_std) {
    return(l%>%
             mutate(res = ifelse(is.na(std_name),
                                 syn,
                                 std_name))%>%
             group_by(id)%>%
             summarise(res = str_flatten(res%>%unique, collapse = '|'))%>%
             pull(res))
  }
  return(l%>%
           mutate(res = ifelse(is.na(sys_name),
                               syn,
                               sys_name))%>%
           group_by(id)%>%
           summarise(res = str_flatten(res%>%unique, collapse = '|'))%>%
           pull(res))
}

tighten <- function(x){
  x%>%
    mutate(std_name = convert_names(Confirmed_deletion, to_std = T))%>%
    transmute(sys_name = Confirmed_deletion, std_name, count)
}

"%+%" <- function(...){
  paste0(...)
}

workingDIR <- "C:/Users/zokmi/Desktop/study/coursework/mtx"
setwd(workingDIR)

seqdirs<-list.dirs(path='../'%+%folder_with_fastq,recursive = FALSE)%+%'/artem'
seqfiles<-unlist(map(seqdirs, ~paste(.x,dir(.x)%>%str_subset("output_count.csv$"),sep='/')))

tablenames <- rep(c('d','g'),each=4)%+%rep(1:2,2,each = 2)%+%rep(1:2,4)


read_files <- function(type,abrev,join=FALSE){
  seqdirs<-get('seqdirs',envir = .GlobalEnv)
  tablenames<-get('tablenames',envir = .GlobalEnv)
  
  seqfiles<-unlist(map(seqdirs, ~paste(.x,dir(.x)%>%str_subset(type%+%".csv$"),sep='/')))
  
  chiffre<-tibble(tablename=tablenames,seqfile=seqfiles)
  
  if(join){
    map2(tablenames%+%abrev,seqfiles,
         function(x,y){
           assign(x,
                  read_csv(y,show_col_types = FALSE)%>%
                    tighten%>%
                    mutate(exp=x),
                  envir=.GlobalEnv)
         }
    )
  }
  else{map2(tablenames%+%abrev,seqfiles,
            function(x,y){
              assign(x,
                     read_csv(y,show_col_types = FALSE),
                     envir=.GlobalEnv)
            }
  )
  }
  return(chiffre)
}

chif<-read_files("output_count",'',TRUE)
read_files("output_dm_count",'dm',TRUE)
#read_files("mixaled",'m',TRUE)
#read_files("mixaled_dm",'mdm',TRUE)
read_files("blasted",'b',TRUE)
read_files("blasted_dm",'bdm',TRUE)

tbs<-c(tablenames,
       tablenames %+% 'dm',
       tablenames %+% 'bdm',
       tablenames %+% 'b')

counts <- purrr::reduce(map(tbs,get),full_join)%>%
  pivot_wider(names_from = exp,values_from = count)%>%
  mutate(across(where(is.numeric),
                ~ifelse(is.na(.x),0,.x)))%>%
  mutate(name=ifelse(is.na(std_name),sys_name,std_name))
colnames(counts)

jcounts <- counts%>%dplyr::select(name, where(is.numeric))

jcounts <- jcounts %>% pivot_longer(- name, 
                                    values_to = "count", 
                                    names_to = "Conditions") %>%
  mutate(Experiments = str_sub(Conditions, 1, 2)) %>%
  group_by(name, Experiments) %>%
  summarise(count = sum(count)) %>% 
  ungroup()%>%
  pivot_wider(names_from = "Experiments", values_from= "count")

write_csv(jcounts,seqdirs[1]%+%'/../../jcounts_exp.csv')


library('GGally')

ggcorr(jcounts,
       nbreaks = 6,
       label = TRUE,
       label_size = 4,
       color = "grey50",
       size=2)

#reference NPVs

#npvise<-function(col){
#  col<-ifelse(col!=0,log(col),NA)
#  dens<-density(col,bw=0.25,kernel='gaussian',na.rm=T)
#  mode<-dens$x[which.max(dens$y)]
#  dif<-col-mode
#  sq<-dif*dif
#  NPV<-ifelse(is.na(col),NA,dif/sqrt(1/nrow(na.omit(jcounts))*sum(sq,na.rm=T)))
#  return(NPV)
#}

jcounts

#Set working directory
setwd("C:/Users/zokmi/Desktop/study/coursework/exp")
jcounts <- read_csv("jcounts_exp.csv")


jc <- jcounts%>%
  filter(if_any(where(is.numeric), ~.x!=0))
n = jc%>%dplyr::select(!name)%>%colSums%>%mean
jc<-jc%>%mutate(across(!name, ~ .x / sum(.x) * n))

jc%>%
  filter(d1>0)%>%
  mutate(across(!name, ~ .x + 1))%>%
  mutate(across(!name, ~ .x / d1))%>%
  mutate(across(!name, ~ log2(.x))) -> to_go



to_go$g1%>%hist(breaks = 1000)


GO_norm<-tibble(pvalue = numeric(), Description = character(), condition = character(), overrep = logical(), geneID = character())


for (i in c('g1_d1', 'g2_d2', 'd2_d1', 'g2_g1')){
  condition <- i
  t=str_split_1(condition,'_')[1]
  u=str_split_1(condition,'_')[2]
  to_go %>%
    dplyr::select(name, starts_with(t) | starts_with(u)) %>%
    mutate(log2FoldChange = !!sym(t) - !!sym(u)) %>% 
    mutate(gene = name) %>% 
    dplyr::select(gene, log2FoldChange) -> res
  
  cut<-res%>%filter(!str_detect(gene,'\\|'), log2FoldChange > 0)
  
  gene_list <- cut$log2FoldChange
  names(gene_list) <- cut$gene%>%yeast_names_convert
  gene_list <- sort(gene_list, decreasing = T)
  
  
  norm_gsea <- clusterProfiler::gseGO(
    geneList = gene_list,
    OrgDb = org.Sc.sgd.db,
    keyType = "ENSEMBL",
    ont = "BP")
  
  
  if ( "core_enrichment" %in% names(norm_gsea@result)){
    gene_clust <- norm_gsea@result$core_enrichment%>%
      map(~yeast_names_convert(str_split(., '\\/')%>%
                                 unlist(), to_std = T))%>%
      map(~str_flatten(., collapse = '/'))%>% as.character()
    GO_norm <- norm_gsea@result%>%
      dplyr::select(Description,pvalue)%>%
      mutate(condition = condition, overrep = T, geneID = gene_clust)%>%
      add_row(GO_norm,.)
  }
  
  cut<-res%>%filter(!str_detect(gene,'\\|'), log2FoldChange < 0)
  
  gene_list <- -cut$log2FoldChange
  names(gene_list) <- cut$gene%>%yeast_names_convert
  gene_list <- sort(gene_list, decreasing = T)
  norm_gsea <- clusterProfiler::gseGO(
    geneList = gene_list,
    OrgDb = org.Sc.sgd.db,
    keyType = "ENSEMBL",
    ont = "BP")
  
  
  if ( "core_enrichment" %in% names(norm_gsea@result)){
    gene_clust <- norm_gsea@result$core_enrichment%>%
      map(~yeast_names_convert(str_split(., '\\/')%>%
                                 unlist(), to_std = T))%>%
      map(~str_flatten(., collapse = '/'))%>% as.character()
    GO_norm <- norm_gsea@result%>%
      dplyr::select(Description,pvalue)%>%
      mutate(condition = condition, overrep = F, geneID = gene_clust)%>%
      add_row(GO_norm,.)
  }
  
}

GO_norm%>%mutate(overrep = factor(overrep, levels = c(T,F)))%>%arrange(condition, overrep, pvalue)%>%
write_tsv('./GSEA_exp.txt')

GSEA_exp_to_split <- read_tsv('./GSEA_exp.txt')
files <- GSEA_exp_to_split$condition %>% unique
map(files, ~ write_tsv(GSEA_exp_to_split %>% filter(str_detect(condition, .x)),'./GSE_exp_' %+% .x %+% '.txt'))


#for gprofiler2

GO_report<-tibble(p_value = numeric(), term_name = character(),  condition = character(), intersection = character(), overrep = logical())

for (i in c('g1_d1', 'g2_d2', 'd2_d1', 'g2_g1')){
  condition <- i
  t=str_split_1(condition,'_')[1]
  u=str_split_1(condition,'_')[2]
  to_go %>%
    dplyr::select(name, starts_with(t) | starts_with(u)) %>%
    mutate(log2FoldChange = !!sym(t) - !!sym(u)) %>% 
    mutate(gene = name) %>% 
    dplyr::select(gene, log2FoldChange) -> res

res%>%cut_significant(log2FoldChange, F)->cut

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
              domain_scope = 'custom_annotated', evcodes = TRUE)#можно убрать аннотейтед изъ скопа, и помѣнять fdr на gSCS
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

geneList <- cut$gene%>%convert_names()

GOres <- gost(geneList, organism = 'scerevisiae', ordered_query = T, measure_underrepresentation = F, user_threshold = 0.05, significant = F, correction_method = 'fdr', custom_bg = bgList, sources = c('GO:BP', 'KEGG'), domain_scope = 'custom_annotated', evcodes = TRUE)#можно убрать аннотейтед изъ скопа, и помѣнять fdr на gSCS
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
