## Report Code
library(magrittr)
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)


## Summary of copies per genome
setwd("/Users/nolson/multicopy_taxa_assignment")
## note removing rows with gi value as NA - due to carriage returns????
cluster_df <- read_tsv("refseq_clusters/cluster_genome_count.tsv") %>% filter(!is.na(gi))
copy_count <- cluster_df %>% separate(cluster, c("clus","threshold","id"),  "_",remove = FALSE) %>% filter(threshold == "1.00") %>% group_by(gi) %>% summarize(copies = n())

ggplot(copy_count) + geom_bar(aes(x = copies)) + 
    theme_bw() + labs(y = "Number of Genomes", x = "16S rRNA Gene Copies")

## Clustering Summary Table
### breakdown of number of clusters
### cluster_genome_count
###     cluster - unique cluster id
###     seqs - sequence ids
###     gi - genome id
###     count - number of sequences per genome in the cluster, sequences from the same genome will have the same value

cluster_df %<>% separate(cluster, c("clus","threshold","id"), 
                         "_",remove = FALSE) %>%
    select(-clus)

## Number of sequences per cluster by threshold - cluster size
cluster_df_size <- cluster_df %>% group_by(threshold, id) %>% summarise(size = n()) 


## Number of sequence copies per cluster
cluster_copy_number <- cluster_df %>% group_by(threshold, id, gi, count) %>% summarise(size = n()) 
    # check size should equal count
with(cluster_copy_number, sum(!(count == size)))
## checks out
cluster_copy_number %<>% select(-size) 
cluster_copy_number %>% 
    group_by(threshold, id) %>% 
    summarize(number_sd =sd(count)) %>%
    mutate(mixed_count = ifelse(number_sd == 0 | is.na(number_sd), "Uniform", "Mix")) %>% 
ggplot() + geom_bar(aes(x = threshold, fill = mixed_count), position ="fill") + theme_bw() + labs(fill = "Copy Number", x = "Clustering Threshold", y = "Proportion of Clusters")

cluster_copy_number %>% 
    group_by(threshold, id) %>% 
    summarize(number_sd =sd(count)) %>%
    filter(number_sd > 0 | is.na())
ggplot() + geom_jitter(aes(x = threshold,y = number_sd)) + theme_bw() + labs(x = "Clustering Threshold", y = "Standard Deviation Copies per Genome")

## number of genomes assigned to multiple clusters
cluster_df_ambig <- cluster_copy_number %>% group_by(threshold, gi) %>% 
    summarise(n_assigned = n()) %>% 
    mutate(ambig = ifelse(n_assigned > 1, "multiple","single"))
cluster_df_ambig %>% ggplot() +
    geom_bar(aes(x = threshold, fill = ambig)) + theme_bw() + 
    labs(x = "Clustering Threshold", y = "Number of Genomes")

cluster_df_ambig %>% filter(n_assigned > 1) %>% 
    ggplot() + geom_bar(aes(x = threshold, fill = factor(n_assigned))) + theme_bw() + labs(x = "Clustering Threshold", y = "Number of Genomes", fill = "Number of Assigned Clusters")


## taxonomic breakdown
taxaID <- read_csv("clusters_with_taxonomic_annotation.csv") %>% filter(!is.na(gi))
cluster_df_anno <- taxaID %>% 
    separate(cluster, c("clus","threshold","id"),"_") %>% 
    select(-clus, -seq_names) %>% 
    group_by(threshold, id) %>% 
    gather(taxa_level, taxa, -threshold, -id, -seqs)

cluster_df_anno_ambig <- cluster_df_anno %>% group_by(threshold, id, taxa_level, taxa) %>% 
    summarise(n_assigned = n()) %>%
    mutate(ambig = ifelse(n_assigned > 1, "multiple","single"))

## This plot does not make sense
cluster_df_anno_ambig %>% filter(threshold != "1.00", 
                                 !(taxa_level %in% c("gi","taxid"))) %>% 
    ggplot() + geom_bar(aes(x = threshold, fill = ambig)) + 
    facet_wrap(~taxa_level, ncol = 1, scale = "free_y") + theme_bw() +
    labs(x = "Cluster Threshold", y = "Number of Clusters", fill = "Proportion of Taxa Assigned to Multiple Clusters")

## Number of clusters with multiple taxa
cluster_df_anno_count <- cluster_df_anno %>% group_by(threshold, id, taxa_level) %>% 
    summarise(taxa_count = length(unique(taxa))) %>% 
    mutate(ambig = ifelse(taxa_count > 1, "multiple","single"))

cluster_df_anno_count %>% filter(threshold != "1.00", 
                                 !(taxa_level %in% c("gi","taxid"))) %>% 
    ggplot() + geom_point(aes(x = threshold, y = taxa_count)) + 
    facet_wrap(~taxa_level, ncol = 1, scale = "free_y") + theme_bw() +
    labs(x = "Cluster Threshold", y = "Number of Taxa Per Cluster") + scale_y_log10()

cluster_df_anno_count %>% filter(threshold != "1.00", 
                                 !(taxa_level %in% c("gi","taxid"))) %>% 
    ggplot() + geom_bar(aes(x = threshold, fill = ambig), position = "fill") + 
    facet_wrap(~taxa_level, ncol = 1, scale = "free_y") + theme_bw() +
    labs(x = "Cluster Threshold")

thresholds <- cluster_df_size$threshold %>% unique() %>% as.numeric() %>% .[-15]
cluster_df_anno_count %>% group_by(threshold, taxa_level, ambig) %>% summarise(ambig_count = n()) %>% group_by(threshold,taxa_level) %>% mutate(ambig_prop = ambig_count/sum(ambig_count)) %>% 
    filter(threshold != "1.00", 
               !(taxa_level %in% c("gi","taxid")),
               ambig == "multiple") %>% 
    ggplot() + geom_point(aes(x = as.numeric(threshold), y = ambig_prop, color = taxa_level)) +
    geom_line(aes(x = as.numeric(threshold), y = ambig_prop, color = taxa_level)) +
    scale_x_continuous(breaks = thresholds, 
                       labels = thresholds) +
    theme_bw() + labs(x = "Cluster Threshold", y = "Proportion of Clusters", color = "Taxonomic Level") +
    geom_hline(aes(yintercept = 0.05), linetype = 2, color = "grey40")