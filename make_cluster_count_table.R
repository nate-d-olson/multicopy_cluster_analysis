## Create copy number analysis files
## For each cluster threshold
## format tsv file
## output column description
## gi - source sequence gi
## cluster - cluster label
## count - number of seqs per cluster
library(magrittr)
library(stringr)
library(dplyr)
library(purrr)
library(readr)

## read cluster input file
read_cluster <- function (clust_file, threshold){
  cluster_set <- read_lines(clust_file, skip = 1) %>% 
      as.list() %>% 
      flatmap(str_split,"\t")
  
  names(cluster_set) <- paste0("cluster_",threshold,"_",1:length(cluster_set))
  cluster_names <- cluster_set  %>% map(length)  %>% 
      unlist()  %>% rep(names(cluster_set), .)
  
  cluster_seqs <- cluster_set  %>% as_vector(.type = "character")
  cluster_df <- data_frame(cluster = cluster_names, seqs = cluster_seqs)
}

## parse seq names
extract_gi_seq_name <- function(seq_name){
    seq_name %>% 
        str_extract('gi\\|.*\\|ref') %>% 
        str_replace_all("gi\\|","") %>% str_replace_all("\\|ref","")
}

## format table
format_cluster_table <- function(cluster_df){
    cluster_df %>% 
        mutate(gi = extract_gi_seq_name(seqs)) %>% # extracting gi
        group_by(cluster, gi) %>% mutate(count = n()) # cluster count per genome
}

setwd("/Users/nolson/multicopy_taxa_assignment/refseq_clusters/")
cluster_files <- list.files(".","*.cluster")

cluster_df <- data_frame()
for(i in cluster_files){
    threshold <- str_extract(i, "\\d\\.\\d\\d")
    cluster_df <- read_cluster(i, threshold) %>% 
        format_cluster_table() %>% 
        bind_rows(cluster_df, .)
}

write_tsv(cluster_df, "cluster_genome_count.tsv")
