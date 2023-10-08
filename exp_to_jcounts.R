library(tidyverse)
library(clusterProfiler) #for GSEA
library('org.Sc.sgd.db') #database for Gene Onthology
library(gprofiler2) #for Over Representation Analysis

folder_with_fastq<-'exp'
#table made with https://www.yeastgenome.org/
#contains systematic and some standard gene names.
to <- read_tsv("../tables/tight_orfs.tsv")%>%
  transmute(sys_name = `Gene > Systematic Name`, std_name = `Gene > Standard Name`)

#function that adds standard gene names to count table and reducts column number
to_join <- function(x){
  y <- x%>%
    transmute(sys_name = Confirmed_deletion, count)%>%
    left_join(to, by = c("sys_name"))%>%
    group_by(sys_name, std_name)%>%
    summarise(count = sum(count))%>%
    ungroup()
  y$count[is.na(y$count)] <- 0
  return(y)
}

#table with other synonyms of gene names. https://www.yeastgenome.org/
yeast_genes <- read_tsv('../tables/yeast_gene_names.txt')
#function converts all gene names either to systematic or standard with to_std=T option. 
convert_names <- function(vals, to_std = F){
  if (length(vals) == 0){return(vals)}
  l <- tibble(syn = vals, id = (1:length(vals)))%>%
   #gene names, combined with '|' are converted separately and returned back joined
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

#to_join analogue
tighten <- function(x){
  x%>%
    mutate(std_name = convert_names(Confirmed_deletion, to_std = T))%>%
    transmute(sys_name = Confirmed_deletion, std_name, count)
}

#three following functions are presented as example of "cut_significant" function in ORA (look below)
llk <- function(params, data) {
  mu <- params[1]
  sigma <- params[2]
  
  # Вычисление логарифма правдоподобия
  ll <- sum(dnorm(data, mean = mu, sd = sigma, log = TRUE))
  
  return(-ll)
}
opt_llh<-function(x){
  #m<-locmodes(x,lowsup=min(x),uppsup=max(x))$locations
  #m<-log(m)
  m<-0
  optim(par=c(m,1), llk, data=x)$par
  
}
cut_significant <- function(df, col, log = T){
  vals <- df %>% pull({{col}})
  if (log) {vals <- log(vals)}
  mu_sigma <- opt_llh(vals)
  quantile <- qnorm(0.95, mean = mu_sigma[1], sd = mu_sigma[2])
  df%>%
    filter(vals > quantile)
}

#paste0 allias
"%+%" <- function(...){
  paste0(...)
}

workingDIR <- "C:/path/wd"
setwd(workingDIR)

#folder_with_fastq has folders for each fastq file with count data table in inner output_folder folder.
seqdirs <- list.dirs(path = '../' %+% folder_with_fastq, recursive = FALSE) %+% '/output_folder' #e.g. ../folder_with_fastq/fastq_name/output_folder
seqfiles <- unlist(map(seqdirs, ~paste(.x, dir(.x) %>% str_subset("output_count.csv$"), sep='/'))) # ../folder_with_fastq/fastq_name/output_folder/XYZ_output_count.csv

#variables for the tables
tablenames <- rep(c('d','g'),each=4)%+%rep(1:2,2,each = 2)%+%rep(1:2,4)

#function reads files in 'seqdirs' directories with suffix defined by type, assigns them to variables named by combination of 'tablenames' and abrev. join = T applies to_join
#also the name of a variable is added to 'exp' column(if join=T) and 'chiffre' table is returned containing variables and filenames correspondance.
read_files <- function(type, abrev, join = FALSE){
  seqdirs <- get('seqdirs', envir = .GlobalEnv)
  tablenames <- get('tablenames', envir = .GlobalEnv)
  
  seqfiles <- unlist(map(seqdirs, ~paste(.x, dir(.x)%>%str_subset(type %+% ".csv$"), sep = '/')))
  
  chiffre<-tibble(tablename=tablenames,seqfile=seqfiles)
  
  if(join){
    map2(tablenames %+% abrev, seqfiles,
         function(x, y){
           assign(x,
                  read_csv(y, show_col_types = FALSE)%>%
                    tighten%>%
                    mutate(exp = x),
                  envir = .GlobalEnv)
         }
    )
  }
  else{map2(tablenames %+% abrev, seqfiles,
            function(x, y){
              assign(x,
                     read_csv(y, show_col_types = FALSE),
                     envir = .GlobalEnv)
            }
  )
  }
  return(chiffre)
}

# reading tables with count data
chif <- read_files("output_count", '', TRUE) # counts from reads with strict barcode and adapter math
read_files("output_dm_count", 'dm', TRUE) # counts from reads with strict barcode and not-strict adapter math
read_files("blasted", 'b', TRUE) # counts from reads with not-strict barcode and strict adapter math
read_files("blasted_dm", 'bdm', TRUE) # counts from reads with not-strict barcode and not-strict adapter math

tbs<-c(tablenames,
       tablenames %+% 'dm',
       tablenames %+% 'bdm',
       tablenames %+% 'b')
#making jcounts table
counts <- purrr::reduce(map(tbs, get), full_join)%>%
  pivot_wider(names_from = exp, values_from = count)%>%
  mutate(across(where(is.numeric),
                ~ifelse(is.na(.x), 0, .x)))%>%
  mutate(name = ifelse(is.na(std_name), sys_name, std_name))
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

#correlations
library('GGally')

ggcorr(jcounts,
       nbreaks = 6,
       label = TRUE,
       label_size = 4,
       color = "grey50",
       size=2)

