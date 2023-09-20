#unannotated part of the script similar to exp_to_GO, but for another data

setwd("C:/path/ngs")
jcounts <- read_csv("jcounts_ngs.csv")


jc <- jcounts%>%
  dplyr::select(name,'pet' %+% 1:3, 'rho1')%>%
  filter(if_any(where(is.numeric), ~.x!=0))
n = jc%>%dplyr::select(!name)%>%colSums%>%mean
jc<-jc%>%mutate(across(!name, ~ .x / sum(.x) * n))

jc%>%
  mutate(across(!name, ~ .x + 1))%>%
  mutate(across(!name, ~ .x / rho1))%>%
  dplyr::select( - rho1)%>%
  pivot_longer(cols = !name, names_to = 'exp')%>%
  mutate(cond = str_sub(exp, 1, 3))%>%
  group_by(name, cond)%>%
  mutate(val = mean(value, na.rm = T))%>%
  ungroup%>%
  dplyr::select(name, cond, val)%>%
  distinct%>%
  pivot_wider(names_from = 'cond', values_from = 'val')%>%
  mutate(across(!name, ~log2(.x)))%>%
  mutate(rho=0)->to_go


to_go$pet%>%hist(breaks = 50)


to_go %>%
    mutate(log2FoldChange = pet - rho) %>% 
    mutate(gene = name) %>% 
    dplyr::select(gene, log2FoldChange) -> res
  
GO_report<-tibble(p_value = numeric(), term_name = character(), intersection = character(), overrep = logical())
  
res %>% filter(log2FoldChange > 0) -> cut
    
geneList <- cut$gene%>%convert_names()
bgList <- res$gene%>%convert_names()
   
GOres <- gost(geneList, organism = 'scerevisiae', ordered_query = T, measure_underrepresentation = F, user_threshold = 0.05, significant = F, correction_method = 'fdr', custom_bg = bgList, sources = c('GO:BP', 'KEGG'), domain_scope = 'custom_annotated', evcodes = TRUE)#можно убрать аннотейтед изъ скопа, и помѣнять fdr на gSCS
GO_report<-GOres$res%>%dplyr::select(p_value, term_name, intersection)%>%as_tibble()%>%mutate(intersection = str_replace_all(intersection, ',', '|')%>%convert_names(to_std = T))%>%arrange(p_value)%>%mutate(overrep = T)%>%add_row(GO_report,.)
    
    
res%>%filter(log2FoldChange<0)->cut

geneList <- cut$gene%>%convert_names()

GOres <- gost(geneList, organism = 'scerevisiae', ordered_query = T, measure_underrepresentation = F, user_threshold = 0.05, significant = F, correction_method = 'fdr', custom_bg = bgList, sources = c('GO:BP', 'KEGG'), domain_scope = 'custom_annotated', evcodes = TRUE)#можно убрать аннотейтед изъ скопа, и помѣнять fdr на gSCS
GO_report<-GOres$res%>%dplyr::select(p_value, term_name, intersection)%>%as_tibble()%>%mutate(intersection = str_replace_all(intersection, ',', '|')%>%convert_names(to_std = T))%>%arrange(p_value)%>%mutate(overrep = F)%>%add_row(GO_report,.)

write_tsv(GO_report, './ngs_GO_??.txt')
