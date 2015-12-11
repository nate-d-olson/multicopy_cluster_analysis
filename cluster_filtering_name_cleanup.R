################################################################################
############## Matching tree id names with cluster sequence names ##############
##
## Duplicate sequences (sequences with same names) were removed prior to multiple
## sequence alignment, and replicte (sequences with identical alignments) were
## removed prior to tree construction.  The following creates a new text file
## with sequence names, new sequence names, and a representative name for
## sequences not in the phylogenetic tree. Duplicate sequences are defined in the
## reduced column (0) and dereplicated sequences (1).
##
################################################################################
library(Biostrings)
library(phylotools)
library(magrittr)
library(dplyr)
library(readr)
library(tidyr)

# Load alignment file prior to name change and get sequence names
old_names_aln <- readDNAMultipleAlignment("infernal_filtered/cmd_line_infernal_1.1_pfiltered.fasta")
old_names <- rownames(seq_aln)

# Load seq file after name change get sequence names

new_names_aln <- read.phylip("infernal_filtered/cmd_line_infernal_1.1.pfiltered.phylip")
new_names <- new_names_aln %>% phy2dat() %>% .$seqNam

# Create data frame with original and new names
name_key <- data_frame(old_names, new_names)

# Load reduced file and annotate dataframe with removed seqs
reduced_names_aln <- read.phylip("infernal_filtered/cmd_line_infernal_1.1.pfiltered.phylip.reduced")
reduced_names <- reduced_names_aln %>% phy2dat() %>% .$col1

name_key_reduced <- reduced_names %>% 
    {data_frame(new_names = ., reduced = 1)} %>% 
    right_join(name_key) %>% 
    mutate(reduced = ifelse(is.na(reduced), 0, 1))

# Load original seq file and get sequence name
cluster_df <- read_tsv("refseq_clusters/cluster_genome_count.tsv") %>% filter(!is.na(gi)) %>% separate(cluster, c("clus","threshold","id"),  "_",remove = FALSE) %>% select(-clus) %>% filter(threshold == "1.00")
cluster_df2 <- cluster_df %>% 
    left_join(name_key_reduced, by = c("seqs" = "old_names")) %>% 
    mutate(dedup = ifelse(is.na(new_names),1,0))

## check to make sure dedup is correctly annodated
cluster_df2  %>% group_by(reduced, dedup)  %>% summarise(count = n())



## Getting representative IDs for 1.00 thresholds, use to annotate seqs no in tree with tree ids
deduped_ids <- cluster_df2 %>% filter(dedup == 0)  %>% 
    group_by(cluster) %>% sample_n(1) %>% 
    rename(rep_name = new_names) %>% select(cluster, rep_name)

## Note sure how to best do this ...
## Maybe use annotated clusters with species taxid
cluster_df3 <- cluster_df2 %>% left_join(deduped_ids) %>% select(seqs, new_names, reduced, dedup, rep_name)

write_tsv(cluster_df3, "name_key.tsv")
