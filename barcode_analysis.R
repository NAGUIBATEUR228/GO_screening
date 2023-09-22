#https://www.bioconductor.org/
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)
library(BioCircos)
library(R.utils)
library(Rbowtie2)
library(Rsamtools)
library(Biostrings)
library(ShortRead)

#folder name on upper level
folder<-'ngs'
#files in inner folders
typefq<-'.fastq'
#remove to shorten file names
excess<-''

#Set working directory
workingDIR <- "C:/Users/zokmi/Desktop/study/coursework/mtx"
setwd(workingDIR)

# reverse_complement_string
rev_compl_string <- function(x) {
  as.character(reverseComplement(DNAStringSet(x)))
}

#read file with reference reads file generated from http://csbio.cs.umn.edu/BEAN-counter/gene_barcode_files/ probable original YKOC table
ref_nm<-read_tsv("../ref_nm.txt")%>%
  transmute(Confirmed_deletion=Original_ORF,
            UPTAG_seqs=Barcode,
            UPTAG_notes='YKOC')

#read file with Reference reads file generated from Puddu et al. 2019
reference <- read_tsv("../puddu.txt")%>%
  select("Confirmed_deletion","UPTAG_seqs","UPTAG_seqs1",'UPTAG_notes')%>%
  filter(!is.na(UPTAG_seqs)|!is.na(UPTAG_seqs1))%>%
  unique()%>%
  pivot_longer(starts_with("UPTAG_seqs"),values_to = "UPTAG_seqs")%>%
  select(!'name')%>%
  filter(!is.na(UPTAG_seqs))%>%unique()%>% #filter(if_any(everything(),~str_detect(.x,"\\|")))
  full_join(.,ref_nm)%>%
  pivot_longer(cols=c(Confirmed_deletion,UPTAG_notes))%>%
  pivot_wider(values_fn = ~str_flatten(.x%>%unique,collapse = '|'))%>%
  mutate(UPTAG_notes=ifelse(is.na(UPTAG_notes),'puddu',UPTAG_notes))
  
write_csv(reference,'../reference.txt')

u1 <- DNAString("GATGTCCACGAGGTCTCT")
u1_rc <- reverseComplement(u1)

u2 <- DNAString("CGTACGCTGCAGGTCGAC")
u2_rc <- reverseComplement(u2)

"%+%" <- function(...){
  paste0(...)
}

fastqdirectories<-list.dirs(path='../'%+%folder,recursive=FALSE)

fastqfiles<-unlist(map(fastqdirectories, ~paste(.x,dir(.x)%>%str_subset(typefq%+%'$'),sep='/')))

for (j in fastqfiles){
setwd(workingDIR)
strm <- FastqStreamer(j)

out_table <- tibble(seq = character(), qual=numeric(), count = integer())

repeat {
#if(FALSE){
  fq <- yield(strm)
  if (length(fq) == 0) {break}
  quals<-quality(fq)
  t1 <- add_column(sread(fq)%>%as.data.frame()%>%as_tibble(),as(quals,"matrix")%>%rowMeans())## consider increase/decrease n here! in tables function
  ## n limits the amplicons diversity
  #n must be large as possible. depends on read number in fq
  colnames(t1)<-c('seq_new','qual')
  attach_this_table <- t1%>%
    group_by(seq_new)%>%
    summarise(qual_new=sum(qual)/n(),
              count_new=n())%>%
    ungroup()
  
  out_table <- out_table %>% 
    full_join(attach_this_table, by = c("seq" = "seq_new"))%>%
    transmute(seq=seq,
              count=replace_na(count,0)+replace_na(count_new,0),
              qual=(replace_na(count,0)*replace_na(qual,0)+replace_na(count_new,0)*replace_na(qual_new,0))/(replace_na(count,0)+replace_na(count_new,0)))
}

#system('python ../fastqstreamer.py '%+%j)
raw_reads_count <- out_table%>%arrange(desc(count))
#mkdirs(j%+%'/../artem/')
dir.create(j%+%'/../artem/',showWarnings = FALSE)
setwd(j%+%'/../artem/')
path_excess<-'../'%+%folder%+%'/'%+%str_match(j,folder%+%'/(.*?)/')[,2]
i<-str_remove(j,path_excess)

write_csv(raw_reads_count, 
                paste(na.omit(str_match(i, '/(.*?)'%+%typefq)[,2])[1], 
                                 "raw_reads_count.csv", 
                      sep = "_"))

barcoded_reads <- raw_reads_count%>% 
  filter(
    str_detect(seq, u1%+%"\\s*(.*?)\\s*"%+%u2)|
      str_detect(seq, u2_rc%+%"\\s*(.*?)\\s*"%+%u1_rc)) %>%
  mutate(barcode = str_match(seq, u1%+%"\\s*(.*?)\\s*"%+%u2)[,2]) %>%
  mutate(barcode = ifelse(
    is.na(barcode),
    str_match(seq,u2_rc%+%"\\s*(.*?)\\s*"%+%u1_rc)[,2],
    barcode))%>%
  filter(!is.na(barcode))%>%#if u1u2? 
  mutate(barcode = ifelse(str_detect(seq, "GATGTCCACGAGGTCTCT"),
                          barcode, 
                          rev_compl_string(barcode)))

write_csv(raw_reads_count%>% 
            filter(!(
              str_detect(seq, u1%+%"\\s*(.*?)\\s*"%+%u2)|str_detect(seq, u2_rc%+%"\\s*(.*?)\\s*"%+%u1_rc))), 
          paste(na.omit(str_match(i,'/(.*?)'%+%typefq)[,2])[1], 
                  "not_barcoded_raw_qual_count.csv", 
                  sep = "_"))

final_result <- barcoded_reads %>% 
  mutate(total_qual=count*qual)%>%
  group_by(barcode) %>%
  summarise(n = n(), count = sum(count),qual=sum(total_qual)/sum(count)) %>%
  ungroup()%>%
  arrange(desc(count))#column n means number of unique reads with barcode

#table with information about quality of barcodes
write_csv(final_result, 
          paste(na.omit(str_match(i,'/(.*?)'%+%typefq)[,2])[1], 
                  "barcoded_qual_count.csv", 
                  sep = "_"))


output <- reference %>% 
  left_join(.,final_result, by = c("UPTAG_seqs" = "barcode"))%>% 
  filter(!is.na(count))%>%#some barcodes from final result are not matched
  group_by(Confirmed_deletion)%>%
  summarise(barcode=str_flatten(unlist(str_split(UPTAG_seqs,'\\|'))%>%unique,collapse = '|'),n=sum(n),count=sum(count),notes=str_flatten(unlist(str_split(UPTAG_notes,'\\|'))%>%unique,collapse = '|'))%>%
  ungroup()%>%
  arrange(desc(count))

write_csv(final_result%>%
  filter(!(barcode %in% reference$UPTAG_seqs))%>%
    arrange(desc(count)),
  paste(na.omit(str_match(i,'/(.*?)'%+%typefq)[,2])[1], 
          "not_matched.csv", sep = "_"))

hist<-output %>% 
  ggplot() +
  aes(x = count)+
  geom_histogram(
    binwidth = 10,
    fill='#F06C98',
    color='black')+
  theme_bw()+
  labs(x = "Количество чтений баркода", y="Количество баркодов с таким числом чтений")+
  scale_y_log10()

ggsave(paste(na.omit(str_match(i,'/(.*?)'%+%typefq)[,2])[1], 
      "hist.png", sep = "_"),
hist,
width = 1280,
height = 720,
units = "px",scale=2)

write_csv(output,
      paste(na.omit(str_match(i,'/(.*?)'%+%typefq)[,2])[1], 
            "output_count.csv", sep = "_"))

}
