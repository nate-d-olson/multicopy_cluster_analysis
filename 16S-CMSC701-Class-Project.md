### Workflow
* get bacterial genome database
* extract 16S using mothur and full 16S primer seqs
* phylogenetic analysis to identify paralogs and orthologs??
* annotate sequences with orthology and paralogy information...
    - how to define/ annotate clusters
    - table with annotation information
        + for each genome - number of copies and cluster assignment
* generate simulated data set for method development
    - Error free data
        + select random subset of genomes
        + extract sequences using different 16S primers
            * full seq
            * commonly used regions
    - Simulated Sequence data
        + same set of random subset of genomes
        + using read simulator same primer sets
    - Real data
        + mock community with seq data - parallel simulated and error free data
* cluster
    - dnaclust datasets with whole genome copy number database
        + different cutoff levels
    - generate count table, id paralog clusters
* taxonomic assignment
    - paralog informed assignment
    - maximization-expectation????
        + cluster reassignment
        + construct new taxonomy based clusters - presence absence


### things to think about
* selection of genomes with different levels of similarity
    - full 16S clustering using different thresholds to identify ambiguous clusters
    - use taxa in ambiguous clusters to define simulated datasets
* Definition of unrealted taxa in the same cluster
    - proportion of clusters with multiple genomes from the same taxonomic level, looking for outliers ....

### To Do
* review copy number correction tools

### Project report
* taxonomic composition of genomes
* number of sequences per taxa
* copy number distribution
* diversity distribution
* classification ambiguity - database, error free, simulated, mock
    - % of seqs in ambiguous clusters
    - % of seqs in incorrect clusters
    - % of seqs in unambiguous and correct clusters
* classification ambiguity after paralog cluster correction

#### hypothesis
* clusters with sequences from different/ unrelated taxa
* multiple clusters with sequences from the same genome
    - to what constraints are these satisfied

## OTU genome assignment and copy number count
* download genomes - look at id values, taxa id values
* extract 16S using mothur
* link gene copy to genome
* link genome to taxonomy - use sqlite db to get taxonomic levels

## 10/21/2015
* genome seqs from genome_purity/results, all _E. coli_
* code running from within `from within /Users/nolson/multicopy_taxa_assignment/ref_genomes`

```
cp /Volumes/Transcend/genomic_purity/results/simdir/*/*fasta.gz .
gunzip *
```

* using rnammer to extract full length 16S sequences from genomes
* run on _E. coli_ genome dataset using `extract_16S.sh` bash script
```
for fa in ../ref_genomes/*fasta;
do
    prefix=$(echo $fa | sed 's/.fasta//')
    rnammer -S bacterial -m ssu -xml $prefix.xml -gff $prefix.gff -f $prefix.fasta < $fa
done
```


* Code to download reference genomes from GenBank
```
#!/usr/bin/bash


#1 bacteria genomes from refseq

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/summary.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/all.fna.tar.gz 

#2 bacteria/archaea genomes from NCBI WGS. Downloaded: 105896 files, 51G in 52m 4s (16.6 MB/s)

wget -r --no-parent --reject "index.html*" ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria_DRAFT/

#3 download fungi genomes from refseq

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/


#4 protozoa genomes from refseq

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/


#5 virus genomes



#6 download taxonomy files
```

* processing all refseq genomes
```
extract_16S(){
    prefix=$(echo $fa | sed 's/.fna//')
    if [ ! -e "$prefix.fasta" ]
    then
        echo $prefix
        rnammer -S bacterial -m ssu -xml $prefix-16S.xml -gff $prefix-16S.gff -f $prefix-16S.fasta < $fa
    fi
}

N=8
(
for fa in */*fna; do 
   ((i=i%N)); ((i++==0)) && wait
   extract_16S "$fa" & 
done
)
```

* clustering refseq 16S
    - note did not do any filtering on input 16S
    - clustered at otu cutoffs levels used by greengenes
    - used script cluster_refseq16S.sh
```
for i in 1.00 0.99 0.97 0.94 0.91 0.88 0.85 0.82 0.79 0.76 0.73 0.70 0.67 0.64 0.61;
do 
    ../dnaclust_repo_release3/dnaclust -d -l -t 8 -e 999 -i all_refseq_16S.fasta -s $i > all_refseq_16S_$i.cluster
done
```

## 11/4/2015
* worked out methods for getting sequence lineage and annotating with cluster assignments - in rmarkdown file "~/Google\ Drive/Courses/CMSC701/16S_MultiCopy_Project/16S-multicopy-cluster-annotation.Rmd"
* next steps generate summary plots and work out script for processing multiple cluster sets
* start thinking about how to format the data for assignments

## 11/15/2015
* worked on manuscript outline - in google docs
* multiple sequence using qiime align_seqs.py
    - Alignments
        + `align_seqs.py -i refseq_clusters/all_refseq_16S.fasta -a pair_hmm`
            * `rRNA_gi|126207488|ref|NC_009053.1|_278448-280265_DIR+` present in DB twice - manually removed 
        + Removed two additional seqs - forgot to record name - get with git entry  
        + And `rRNA_gi|386853410|ref|NC_017717.1|_446997-448520_DIR-`
            * Running in parallel
                - `parallel_align_seqs_pynast.py -i refseq_clusters/all_refseq_16S.fasta -a pair_hmm -O 4 -o pynast_aligned_parallel &`
                - Did not complete overnight running non-parallel
                    + `align_seqs.py -i remove_replicate_seqs/no_dups.fasta -a pair_hmm &`
        + Duplicates removed using perl one liner
            * `perl -ne 'BEGIN{$/=">";$"=";"}($d,$_)=/(.*?)\n(.+?)>?$/s;push @{$h{lc()}},$d if$_;END{for(keys%h){print">@{$h{$_}}\n$_"}}' all_refseq_16S.fasta > no_dups.fasta`
                - modified from https://www.biostars.org/p/3003/
                - Restarted pynast alignment - `parallel_align_seqs_pynast.py -i remove_replicate_seqs/no_dups.fasta -a pair_hmm -O 8 -o pynast_aligned_parallel &`

        + Infernal
            * downloaded reference alignment from http://rfam.xfam.org/family/SSU_rRNA_bacteria?tab=alignBlock#tabview=tab2
            * `align_seqs.py -i refseq_clusters/all_refseq_16S.fasta -m infernal -t infernal/RF00177.stockholm.txt -o infernal_aligned &`
            * Running from the command line outside of pynast
            * Error with mxsize
            ```
            cmalign --verbose --ifile info.txt --elfile seq_el.txt --sfile score.txt --tfile parsetrees.txt -o ../infernal_aligned/cmd_line_infernal_1.1.fasta --cpu 7 --mxsize 5000 --informat FASTA --dnaout RF00177.cm ../remove_replicate_seqs/no_dups.fasta
            ```
            renamed to cmd_line_infernal1.1.fasta to_cmd_line_infernal_1.1.sto to
            used to convert stokholm format to fasta http://www.hiv.lanl.gov/content/sequence/FORMAT_CONVERSION/form.html
    - generate phylogenetic tree and calculate distances based using methods from Eisen's group
        + " Reference sequences were aligned to the GreenGenes core
 set   [29]   with   PyNAST   [30]   and   masked   with   the   GreenGenes
 lanemask. We constructed a phylogeny of the reference sequences
 using RAxML [31] with a GTR+Gamma model of evolution.
"
        + Incorporating 16S Gene Copy Number Information
 Improves Estimates of Microbial Diversity and
 Abundance
    - use filter alignment - `http://qiime.org/scripts/filter_alignment.html`
        - `filter_alignment.py -i infernal_aligned/cmd_line_infernal_1.1.fasta -o filtered_infernal/`
        - issues with file format ....
    ad make tree in qiime ....

* catadapt - more flexibility in adaptor trimming
* publications on the impact of trimming and quality filtering
* less -S // no wrapping 

Infernal 

### Creating taxonomy database
* using taxtastic python package to install 
    - http://fhcrc.github.io/taxtastic
* installed using `pip install taxtastic`
* command line 
    - `taxit new_database -d ncbi_taxa_db.sqlite`
* Adding gi table
    - downloaded `gi_taxid_nucl.dmp.gz` from `ftp://ftp.ncbi.nih.gov/pub/taxonomy`
    - sqlite script `make_tax_db`
    - command line `sqlite3 ncbi_taxa_db.sqlite < ../make_tax_db`

### generating multiple sequence alignment with infernal
* removed duplicate seqs with - `2015-11-30-16S_copy_find_duplicates.Rmd`
    * manually removed identified replicates
+ running infernal from the command line started at 9:18 PM
    * computer crashed had to restart 9:49 completed 10:37 PM
    * from within `infernal_aligned`
    * `cmalign --verbose --ifile info.txt --elfile seq_el.txt --sfile score.txt --tfile parsetrees.txt -o ../infernal_aligned/16S_infernal.fasta --cpu 7 --mxsize 5000 --informat FASTA --dnaout RF00177.cm ../refseq_genomes/all_refseq_16S.fasta`
+ Filtering alignment started 10:40 PM
    * - `filter_alignment.py -i infernal_aligned/cmd_line_infernal_1.1.fasta -o infernal_filtered/ -r`
    * added option -r to remove outliers as RAXML had error for sequences with missing data
    * using -r did not remove seq manually removed sequences
        - `>rRNA_gi|525896797|ref|NC_021827.1|_2952632-2953716_DIR-`
        - `>rRNA_gi|408675720|ref|NC_018750.1|_3500124-3501493_DIR-`
        - `>rRNA_gi|560141539|ref|NC_022904.1|_2-1276_DIR+`
        - `>rRNA_gi|533108975|ref|NC_022132.1|_1583351-1584486_DIR`
```
raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -T 7 -n infernal_16S_tree -s infernal_filtered/cmd_line_infernal_1.1.pfiltered.phylip
ERROR: Sequence ID.1512 consists entirely of undetermined values which will be treated as missing data
ERROR: Sequence ID.1749 consists entirely of undetermined values which will be treated as missing data
ERROR: Sequence ID.4742 consists entirely of undetermined values which will be treated as missing data
ERROR: Sequence ID.4755 consists entirely of undetermined values which will be treated as missing data
```



+ Tree generating code from eisen study
    * converted to phylip with easel
        - `esl-reformat --rename ID phylip infernal_filtered/cmd_line_infernal_1.1_pfiltered.fasta > infernal_filtered/cmd_line_infernal_1.1.pfiltered.phylip`
    * `raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -T 7 -p 16 -n infernal_16S_tree -s infernal_filtered/cmd_line_infernal_1.1.pfiltered.phylip`
        - started at ~11:45 
        - error messages about identical sequence names modified code 
    * completed after about 3 hours
    * moved files to `infernal_raxml_tree`
+ Next step - generate distance matrix from tree using `cophenetic.phylo` in the R ape package
