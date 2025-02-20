library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(BSgenome)
library(biomaRt)



# get gene name of transcripts
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
chromosome <- c(1:22, "X", "Y")
# biomaRt_list <- listFilters(ensembl)
for (chr in chromosome){
  df <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", "refseq_mrna", "transcript_mane_select"),
              filters = "chromosome_name", values = chr,
              mart = ensembl)
  
  if(chr == chromosome[1]){
    gene_name_of_transcript <- df
  } else {
    gene_name_of_transcript <- bind_rows(gene_name_of_transcript, df)
  }
}
gene_name_of_transcript <- gene_name_of_transcript %>% 
  mutate(gene_name = external_gene_name) %>% 
  mutate(gene_id = ensembl_gene_id) %>% 
  mutate(transcript_id = ensembl_transcript_id) %>% 
  dplyr::select(-external_gene_name, -ensembl_gene_id, -ensembl_transcript_id)



# get all exons and introns from hg38
getBSgenome("BSgenome.Hsapiens.UCSC.hg38", masked=FALSE, load.only=FALSE)

exons_df <- data.frame(exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "tx", use.names = TRUE)) %>% 
  mutate(transcript_id = group_name) %>% 
  #dplyr::select(-group_name) %>% 
  mutate(transcript_id = sub("\\..*", "", transcript_id)) %>% 
  left_join(gene_name_of_transcript, by = "transcript_id")

introns_df <- data.frame(intronsByTranscript(TxDb.Hsapiens.UCSC.hg38.knownGene, use.names = TRUE)) %>% 
  mutate(transcript_id = group_name) %>% 
  #dplyr::select(-group_name) %>% 
  mutate(transcript_id = sub("\\..*", "", transcript_id)) %>% 
  left_join(gene_name_of_transcript, by = "transcript_id")



# export to CSV files
write_csv(exons_df, "exons_in_hg38.csv.gz")
write_csv(introns_df, "introns_in_hg38.csv.gz")
